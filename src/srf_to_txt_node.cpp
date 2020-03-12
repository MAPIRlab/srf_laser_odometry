/* This script prints to file the GT pose and estimations of SRF
 *
 * */
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <mrpt_bridge/pose.h>
#include <mrpt/poses/CPose3D.h>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;


//Globals
bool hasGT, hasSRF_CS, hasSRF_KS, hasSRF_HYBRID, hasCSM, hasPSM;
mrpt::poses::CPose3D GTpose, SRF_CSpose, SRF_KSpose, SRF_HYBRIDpose, CSMpose, PSMpose;
double srf_freq;
std::ofstream myfile;
bool first_estimation = true;



//Callback every time a new GT pose is available
void GTOdomCallBack(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, GTpose);
    hasGT = true;
}


//Callback every time a SRF_CS estimation is published!
void srfOdomCallBack_cs(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, SRF_CSpose);
    hasSRF_CS = true;
}

//Callback every time a SRF_KS estimation is published!
void srfOdomCallBack_ks(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, SRF_KSpose);
    hasSRF_KS = true;
}

//Callback every time a SRF_HYBRID estimation is published!
void srfOdomCallBack_hybrid(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, SRF_HYBRIDpose);
    hasSRF_HYBRID = true;
}

//Callback every time a CSM estimation is published!
void csmOdomCallBack(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, CSMpose);
    hasCSM = true;
}

//Callback every time a SRF_CS estimation is published!
void psmOdomCallBack(const nav_msgs::Odometry::ConstPtr& new_odom)
{
    mrpt_bridge::convert(new_odom->pose.pose, PSMpose);
    hasPSM = true;
}


void writeToTextFile()
{
    if (hasGT && hasSRF_CS && hasSRF_KS && hasSRF_HYBRID && hasCSM && hasPSM)
    {
        // Freq
        myfile << std::fixed << std::setprecision(2) << srf_freq  << " ";
        //myfile << std::fixed << std::setprecision(8) << ros::Time::now().toSec() << " ";

        //GT pose (most actual one - may not be the closest one)
        myfile << std::fixed << std::setprecision(8) << GTpose.x()   << " ";
        myfile << std::fixed << std::setprecision(8) << GTpose.y()   << " ";
        myfile << std::fixed << std::setprecision(8) << GTpose.yaw() << " ";

        //SRF_CS estimation
        myfile << std::fixed << std::setprecision(8) << SRF_CSpose.x()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_CSpose.y()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_CSpose.yaw() << " ";

        //SRF_KS estimation
        myfile << std::fixed << std::setprecision(8) << SRF_KSpose.x()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_KSpose.y()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_KSpose.yaw() << " ";

        //SRF_HYBRID estimation
        myfile << std::fixed << std::setprecision(8) << SRF_HYBRIDpose.x()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_HYBRIDpose.y()   << " ";
        myfile << std::fixed << std::setprecision(8) << SRF_HYBRIDpose.yaw() << " ";

        //CSM estimation
        myfile << std::fixed << std::setprecision(8) << CSMpose.x() << " ";
        myfile << std::fixed << std::setprecision(8) << CSMpose.y() << " ";
        myfile << std::fixed << std::setprecision(8) << CSMpose.yaw() << " ";

        //PSM estimation
        myfile << std::fixed << std::setprecision(8) << PSMpose.x() << " ";
        myfile << std::fixed << std::setprecision(8) << PSMpose.y() << " ";
        myfile << std::fixed << std::setprecision(8) << PSMpose.yaw();

        myfile << "\n";

        //reset all
        hasGT = hasSRF_CS = hasSRF_KS = hasSRF_HYBRID = hasCSM = hasPSM = false;
    }
    else
    {
        //ROS_INFO("GT: %s", hasGT ? "true" : "false");
        //ROS_INFO("CS: %s", hasSRF_CS ? "true" : "false");
        //ROS_INFO("KS: %s", hasSRF_KS ? "true" : "false");
        //ROS_INFO("HY: %s", hasSRF_HYBRID ? "true" : "false");
        //ROS_INFO("CSM: %s", hasCSM ? "true" : "false");
        //ROS_INFO("PSM: %s\n", hasPSM ? "true" : "false");
    }
}





// MAIN
int main(int argc, char** argv)
{
    //Init
    ros::init(argc, argv, "ros2txt_node");
    ros::NodeHandle n;
    ros::NodeHandle pn("~");
    hasGT = hasSRF_CS = hasSRF_KS = hasSRF_HYBRID = hasCSM = hasPSM = false;

    //Params
    std::string GT_odom_topic = "/robot_0/base_pose_ground_truth";
    std::string srf_cs_odom_topic = "/odom_srf_cs";
    std::string srf_ks_odom_topic = "/odom_srf_ks";
    std::string srf_hybrid_odom_topic = "/odom_srf_hybrid";
    std::string csm_odom_topic = "/odom_csm";
    std::string psm_odom_topic = "/odom_psm";

    n.getParam("/srf_laser_odometry/laser_decimation", srf_freq);
    srf_freq = 10/srf_freq;   //Chapuza suponiendo que el laser va a 10Hz

    //Publishers & subscriptors
    ros::Subscriber GTOdom_sub = n.subscribe<nav_msgs::Odometry>(GT_odom_topic,1,&GTOdomCallBack);
    ros::Subscriber srf_sub_cs = n.subscribe<nav_msgs::Odometry>(srf_cs_odom_topic,1,&srfOdomCallBack_cs);
    ros::Subscriber srf_sub_ks = n.subscribe<nav_msgs::Odometry>(srf_ks_odom_topic,1,&srfOdomCallBack_ks);
    ros::Subscriber srf_sub_hybrid = n.subscribe<nav_msgs::Odometry>(srf_hybrid_odom_topic,1,&srfOdomCallBack_hybrid);
    ros::Subscriber csm_sub = n.subscribe<nav_msgs::Odometry>(csm_odom_topic,1,&csmOdomCallBack);
    ros::Subscriber psm_sub = n.subscribe<nav_msgs::Odometry>(psm_odom_topic,1,&psmOdomCallBack);


    //Open file to save data
    char buffer [50];
    std::sprintf(buffer, "/home/jgmonroy/poses_srf_to_txt_%.1fHz.txt", srf_freq);
    myfile.open(buffer);
    myfile << "Freq(Hz) GT_pose(x y thetha) SRF_CS_pose(x y thetha) SRF_KS_pose(x y thetha) SRF_HYBRID_pose(x y thetha) CSM_pose(x y thetha) PSM_pose(x y thetha)\n";


    ROS_INFO("Starting SAVETOTEXT");
    //Loop
    ros::Rate r(100);
    while(ros::ok())
    {
        writeToTextFile();
        ros::spinOnce();
        r.sleep();
    }

    myfile.close();


    return(0);
}



