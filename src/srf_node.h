/** ****************************************************************************************
*  This node presents a fast and precise method to estimate the planar motion of a lidar
*  from consecutive range scans. It is very useful for the estimation of the robot odometry from
*  2D laser range measurements.
*  This module is developed for mobile robots with innacurate or inexistent built-in odometry.
*  It allows the estimation of a precise odometry with low computational cost.
*  For more information, please refer to:
*
*  Planar Odometry from a Radial Laser Scanner. A Range Flow-based Approach. ICRA'16.
*  Available at: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/217
*  Authors: Mariano Jaimez Tarifa, Javier G. Monroy. MAPIR group, University of Malaga, Spain
*
*  Maintainer: Javier G. Monroy
*  MAPIR group: http://mapir.isa.uma.es/
*  Date: January 2016
*  More Info: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/217
*********************************************************************/

#ifndef CLaserOdometry2D_H
#define CLaserOdometry2D_H

#include <ros/ros.h>
#include <stdlib.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <mrpt_bridge/pose.h>

// SRF
#include "laser_odometry_refscans.h"

// MRPT related headers
#include <mrpt/version.h>
#if MRPT_VERSION>=0x130
#	include <mrpt/obs/CObservation2DRangeScan.h>
#   include <mrpt/obs/CObservationOdometry.h>
    using namespace mrpt::obs;
#else
#	include <mrpt/slam/CObservation2DRangeScan.h>
#   include <mrpt/slam/CObservationOdometry.h>
    using namespace mrpt::slam;
#endif
#include <mrpt/poses/CPose3D.h>




class CLaserOdometry2D
{
public:
    //Variables
    std::string laser_scan_topic;
    std::string base_frame_id;
    std::string odom_topic, odom_frame_id;
    std::string init_pose_from_topic;
    std::string operation_mode;
    double laser_min_range, laser_max_range;
    bool publish_tf;
    mrpt::poses::CPose3D robot_pose, robot_oldpose;
    int laser_counter, laser_decimation;

    //Core class of SRF
    SRF_RefS srf_obj;

    //methods
    CLaserOdometry2D();
    ~CLaserOdometry2D();
    bool is_initialized();
    bool scan_available();
    void Init();
    void odometryCalculation();     //Update odometric pose
    void publishPoseFromSRF();     //Publishes the last odometric pose with ROS format

protected:
    ros::NodeHandle n;
    sensor_msgs::LaserScan last_scan;
    bool module_initialized, first_laser_scan, new_scan_available, GT_pose_initialized;
    tf::TransformListener tf_listener;                              //Do not put inside the callback
    tf::TransformBroadcaster odom_broadcaster;
    ros::Time last_odom_time;
    nav_msgs::Odometry initial_robot_pose;

    //Subscriptions & Publishers
    ros::Subscriber laser_sub, initPose_sub;
    ros::Publisher odom_pub, laser_pub;

    //CallBacks
    void LaserCallBack(const sensor_msgs::LaserScan::ConstPtr& new_scan);
    void initPoseCallBack(const nav_msgs::Odometry::ConstPtr& new_initPose);
};

#endif
