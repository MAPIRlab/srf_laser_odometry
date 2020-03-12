/*********************************************************************
*
* Software License Agreement (BSD License)
*
*  Copyright (c)  2015, Örebro University, Sweden
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*  Authors: Mariano Jaimez Tarifa, Javier G. Monroy
*           MAPIR group, University of Malaga, Spain
*  Date: January 2016
*  More Info: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/217
*********************************************************************/

#include "srf_node.h"

using namespace mrpt;
using namespace mrpt::math;
using namespace mrpt::obs;
using namespace mrpt::poses;
using namespace std;
using namespace Eigen;


// --------------------------
// CLaserOdometry2D Wrapper
//---------------------------
CLaserOdometry2D::CLaserOdometry2D()
{
    ROS_INFO("[SRF] Initializing SRF_node...");

    //Read Parameters
    //----------------
    ros::NodeHandle pn("~");
    pn.param<std::string>("laser_scan_topic",laser_scan_topic,"/laser_scan");
    pn.param<bool>("publish_tf", publish_tf, true);
    pn.param<std::string>("base_frame_id", base_frame_id, "/base_link");
    pn.param<std::string>("odom_topic", odom_topic, "/odom");
    pn.param<std::string>("odom_frame_id", odom_frame_id, "/odom");
    pn.param<std::string>("init_pose_from_topic", init_pose_from_topic, "");
    pn.param<int>("laser_decimation",laser_decimation,1);
    pn.param<double>("laser_min_range",laser_min_range,-1.0);
    pn.param<double>("laser_max_range",laser_max_range,-1.0);

    pn.param<std::string>("operation_mode", operation_mode, "HYBRID");  //CS=consecutiveScans, KS=keyScans, HYBRID=threeScansWithKeyScan

    //Publishers and Subscribers
    //--------------------------    
    laser_sub = n.subscribe<sensor_msgs::LaserScan>(laser_scan_topic,1,&CLaserOdometry2D::LaserCallBack,this);
    odom_pub = pn.advertise<nav_msgs::Odometry>(odom_topic, 5);
    laser_pub = pn.advertise<sensor_msgs::LaserScan>("srf_laser_truncated", 5);

    //init pose
    //----------
    if (init_pose_from_topic != "")
    {
        initPose_sub = n.subscribe<nav_msgs::Odometry>(init_pose_from_topic,1,&CLaserOdometry2D::initPoseCallBack,this);
        GT_pose_initialized  =false;
    }
    else
    {
        GT_pose_initialized = true;
        initial_robot_pose.pose.pose.position.x = 0;
        initial_robot_pose.pose.pose.position.y = 0;
        initial_robot_pose.pose.pose.position.z = 0;
        initial_robot_pose.pose.pose.orientation.w = 0;
        initial_robot_pose.pose.pose.orientation.x = 0;
        initial_robot_pose.pose.pose.orientation.y = 0;
        initial_robot_pose.pose.pose.orientation.z = 1;
    }

    //Init variables
    //----------------
    module_initialized = false;
    first_laser_scan = true;
    laser_counter = 0;
}


CLaserOdometry2D::~CLaserOdometry2D()
{
}


bool CLaserOdometry2D::is_initialized()
{
    return module_initialized;
}


bool CLaserOdometry2D::scan_available()
{
    return new_scan_available;
}


void CLaserOdometry2D::Init()
{    
    ROS_INFO("[SRF] Got first Laser Scan .... Configuring node");
    const unsigned int scan_size = last_scan.ranges.size();             // Num of samples (size) of the scan laser    
    const float fov = fabs(last_scan.angle_max - last_scan.angle_min);  // Horizontal Laser's FOV

    //Init core class with Laser specs and corresponding operation_mode
    if (!strcmp(operation_mode.c_str(),"CS"))
        srf_obj.initialize(scan_size, fov, 0);
    else if (!strcmp(operation_mode.c_str(),"KS"))
        srf_obj.initialize(scan_size, fov, 1);
    else if (!strcmp(operation_mode.c_str(),"HYBRID"))
        srf_obj.initialize(scan_size, fov, 2);
    else
    {
        ROS_ERROR("[srf] Operation mode not implemented. Using default HYBRID");
        srf_obj.initialize(scan_size, fov, 2);
    }


    //Set laser pose on the robot (through tF)
    // This allow estimation of the odometry with respect to the robot base reference system.
    mrpt::poses::CPose3D LaserPoseOnTheRobot;
    tf::StampedTransform transform;
    try
    {
        tf_listener.lookupTransform(base_frame_id, last_scan.header.frame_id, ros::Time(0), transform);
    }
    catch (tf::TransformException &ex)
    {
        ROS_ERROR("%s",ex.what());
        ros::Duration(1.0).sleep();
    }

    //TF:transform -> mrpt::CPose3D (see mrpt-ros-bridge)
    const tf::Vector3 &t = transform.getOrigin();
    LaserPoseOnTheRobot.x() = t[0];
    LaserPoseOnTheRobot.y() = t[1];
    LaserPoseOnTheRobot.z() = t[2];
    const tf::Matrix3x3 &basis = transform.getBasis();
    mrpt::math::CMatrixDouble33 R;
    for(int r = 0; r < 3; r++)
        for(int c = 0; c < 3; c++)
            R(r,c) = basis[r][c];
    LaserPoseOnTheRobot.setRotationMatrix(R);
    //LaserPoseOnTheRobot.setYawPitchRoll(LaserPoseOnTheRobot.yaw() - 0.1f, 0.f, 0.f);

    //Robot initial pose
    mrpt::poses::CPose3D robotInitialPose;
    geometry_msgs::Pose _src = initial_robot_pose.pose.pose;
    mrpt_bridge::convert(_src,robotInitialPose);

    //Set the Laser initial pose = Robot_initial_pose + LaserPoseOnTheRobot
    srf_obj.laser_pose = CPose2D(robotInitialPose + LaserPoseOnTheRobot);
    srf_obj.laser_oldpose = CPose2D(robotInitialPose + LaserPoseOnTheRobot);

    module_initialized = true;
    last_odom_time = ros::Time::now();
    ROS_INFO("[SRF] Configuration Done.");
}



void CLaserOdometry2D::publishPoseFromSRF()
{

    ROS_INFO("[SRF] LASERodom = [%f %f %f]",srf_obj.laser_pose.x(),srf_obj.laser_pose.y(),srf_obj.laser_pose.phi());

    // GET ROBOT POSE from LASER POSE
    //--------------------------------
    mrpt::poses::CPose3D LaserPoseOnTheRobot_inv;
    tf::StampedTransform transform;
    try
    {
        tf_listener.lookupTransform(last_scan.header.frame_id, base_frame_id, ros::Time(0), transform);
    }
    catch (tf::TransformException &ex)
    {
        ROS_ERROR("%s",ex.what());
        ros::Duration(1.0).sleep();
    }

    //TF:transform -> mrpt::CPose3D (see mrpt-ros-bridge)
    const tf::Vector3 &t = transform.getOrigin();
    LaserPoseOnTheRobot_inv.x() = t[0];
    LaserPoseOnTheRobot_inv.y() = t[1];
    LaserPoseOnTheRobot_inv.z() = t[2];
    const tf::Matrix3x3 &basis = transform.getBasis();
    mrpt::math::CMatrixDouble33 R;
    for(int r = 0; r < 3; r++)
        for(int c = 0; c < 3; c++)
            R(r,c) = basis[r][c];
    LaserPoseOnTheRobot_inv.setRotationMatrix(R);

    //Compose Transformations
    robot_pose = srf_obj.laser_pose + LaserPoseOnTheRobot_inv;
    ROS_DEBUG("[SRF] BASEodom = [%f %f %f]",robot_pose.x(),robot_pose.y(),robot_pose.yaw());


    // Estimate linear/angular speeds (mandatory for base_local_planner)
    // last_scan -> the last scan received (practically now)
    // last_odom_time -> The time of the previous scan lasser used to estimate the pose
    //-------------------------------------------------------------------------------------
    double time_inc_sec = (last_scan.header.stamp - last_odom_time).toSec();
    last_odom_time = last_scan.header.stamp;
    if (time_inc_sec <=0)
        ROS_WARN("[SRF] Time increment between Odom estimation is: %.6f sec",time_inc_sec);
    else
    {
        double lin_speed_x = srf_obj.kai_loc(0) / time_inc_sec;
        double lin_speed_y = srf_obj.kai_loc(1) / time_inc_sec;
        double ang_speed = srf_obj.kai_loc(2) / time_inc_sec;
        robot_oldpose = robot_pose;

        //first, we'll publish the odometry over tf
        //---------------------------------------
        if (publish_tf)
        {
            geometry_msgs::TransformStamped odom_trans;
            odom_trans.header.stamp = ros::Time::now();
            odom_trans.header.frame_id = odom_frame_id;
            odom_trans.child_frame_id = base_frame_id;
            odom_trans.transform.translation.x = robot_pose.x();
            odom_trans.transform.translation.y = robot_pose.y();
            odom_trans.transform.translation.z = 0.0;
            odom_trans.transform.rotation = tf::createQuaternionMsgFromYaw(robot_pose.yaw());
            //send the transform
            odom_broadcaster.sendTransform(odom_trans);
        }

        //next, we'll publish the odometry message over ROS topic
        //-------------------------------------------------------
        nav_msgs::Odometry odom;
        //odom.header.stamp = ros::Time::now();
        odom.header.stamp = last_scan.header.stamp;
        odom.header.frame_id = odom_frame_id;
        //set the position
        odom.pose.pose.position.x = robot_pose.x();
        odom.pose.pose.position.y = robot_pose.y();
        odom.pose.pose.position.z = 0.0;
        odom.pose.pose.orientation = tf::createQuaternionMsgFromYaw(robot_pose.yaw());
        //set the velocity
        odom.child_frame_id = base_frame_id;
        odom.twist.twist.linear.x = lin_speed_x;    //linear speed
        odom.twist.twist.linear.y = lin_speed_y;
        odom.twist.twist.angular.z = ang_speed;   //angular speed
        //publish the message
        odom_pub.publish(odom);
    }

    //Clear current laser
    new_scan_available = false;
}


//-----------------------------------------------------------------------------------
//                                   CALLBACKS
//-----------------------------------------------------------------------------------

void CLaserOdometry2D::LaserCallBack(const sensor_msgs::LaserScan::ConstPtr& new_scan)
{
    if (GT_pose_initialized)
    {
        //Keep in memory the last received laser_scan
        last_scan = *new_scan;
        laser_counter++;
        //ROS_INFO([SRF] "Laser Counter %u",laser_counter);


        //FOR SIMULATION ONLY
        //-------------------
        //Add noise to scans
        /*
        double noise_mean = 0.0;
        double noise_std = 0.01;
        for (unsigned int i = 0; i<last_scan.ranges.size(); i++)
        {
            //Simualte Gaussian noise
            double x = 0.0;
            for( int j=0; j<12; j++ )
                x += std::rand()/(RAND_MAX+1.0);

            last_scan.ranges[i] += (noise_std*(x-6.0) + noise_mean);
        }
        */


        //Initialize module on first scan
        if (first_laser_scan)
        {
            Init();
            first_laser_scan = false;
            //Load first scan and build pyramid
            for (unsigned int i = 0; i<last_scan.ranges.size(); i++)
                srf_obj.range_wf(i) = last_scan.ranges[i];
            srf_obj.createScanPyramid();

            //Set laser min-max distances (if not set as parameters)
            if (laser_min_range == -1)
                laser_min_range = 0.0;
            if (laser_max_range == -1)
                laser_max_range = 0.99f*last_scan.range_max;
        }
        else
        {
            if (laser_counter % laser_decimation == 0)
            {
                //copy laser scan to internal srf variable
                for (unsigned int i = 0; i<last_scan.ranges.size(); i++)
                {
                    //Check min-max distances, e.g. to avoid including points of the own robot
                    if ( (last_scan.ranges[i] > laser_max_range) || (last_scan.ranges[i] < laser_min_range) )
                    {
                        srf_obj.range_wf(i) = 0.f; //invalid measurement
                        last_scan.ranges[i] = 0.0;
                    }
                    else                    
                        srf_obj.range_wf(i) = last_scan.ranges[i];
                }

                //Publish truncated laser for visualization
                laser_pub.publish(last_scan);

                //Check the average laser measurement
                float sum = 0.f;
                for (unsigned int i = 0; i<last_scan.ranges.size(); i++)
                {
                    sum += srf_obj.range_wf(i);
                }

                ROS_INFO("[SRF] The average laser measurement is %f", sum);

                //Process odometry estimation
                srf_obj.odometryCalculation();
                publishPoseFromSRF();
            }
        }
    }
}


void CLaserOdometry2D::initPoseCallBack(const nav_msgs::Odometry::ConstPtr& new_initPose)
{
    // Initialize robot pose on first GT pose. Else do Nothing!
    // Usefull for comparison with other odometry methods.
    if (!GT_pose_initialized)
    {
        initial_robot_pose = *new_initPose;
        GT_pose_initialized = true;
    }
}


//-----------------------------------------------------------------------------------
//                                   MAIN
//-----------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    ros::init(argc, argv, "SRF_LaserOdom");

    //Wrapper class
    CLaserOdometry2D myLaserOdom;

    //Main Loop
    //----------
    ros::spin();

    return(0);
}
