#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <algorithm>    // std::min

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{

	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// set a reference target speed
double ref_vel = 0.;//49.5; // mph


vector<vector<double>> generate_path(bool too_close, int lane, json j,
  vector<double> map_waypoints_x,  vector<double> map_waypoints_y,
  vector<double> map_waypoints_s){

  double car_x = j["x"];
  double car_y = j["y"];
  double car_s = j["s"];
  double car_yaw = j["yaw"];

  // Previous path data given to the Planner
  auto previous_path_x = j["previous_path_x"];
  auto previous_path_y = j["previous_path_y"];

  int prev_size = previous_path_x.size();

  // could be more efficient to change velocity in path planner.
  if(too_close){
    // ref_vel -= .224;
    ref_vel -= .184;

  }else if(ref_vel < 49.5){
    ref_vel += .224;
  }

  // create a list of widely spaced (x, y) waypoints, evenly spaced at 30m
  // later we will interpolate these waypoints with a spline and fill it
  // with more points that control speed
  vector<double> ptsx;
  vector<double> ptsy;

  // reference x, y, yaw states
  // either we will reference the starting point as where the car is or at the previous paths and point
  double ref_x = car_x;
  double ref_y = car_y;
  double ref_yaw = deg2rad(car_yaw);

  // if previous size is almost empty use the car as starting point
  if(prev_size < 2){
    // use two points that make the path tangent to the car
    double prev_car_x = car_x - cos(ref_yaw);
    double prev_car_y = car_y - sin(ref_yaw);

    ptsx.push_back(prev_car_x);
    ptsx.push_back(car_x);

    ptsy.push_back(prev_car_y);
    ptsy.push_back(car_y);
  }
  else{
    // use the previous path's end point as starting reference
    ref_x = previous_path_x[prev_size-1];
    ref_y = previous_path_y[prev_size-1];

    double ref_x_prev = previous_path_x[prev_size-2];
    double ref_y_prev = previous_path_y[prev_size-2];
    ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

    ptsx.push_back(ref_x_prev);
    ptsx.push_back(ref_x);

    ptsy.push_back(ref_y_prev);
    ptsy.push_back(ref_y);

  }

  // In Frenet coordinates, add evenly 30m points ahead of the starting reference point
  vector<double> next_wp0 = getXY(car_s+50, (2+4*lane), map_waypoints_s,map_waypoints_x,map_waypoints_y);
  vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s,map_waypoints_x,map_waypoints_y);
  vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s,map_waypoints_x,map_waypoints_y);

  ptsx.push_back(next_wp0[0]);
  ptsx.push_back(next_wp1[0]);
  ptsx.push_back(next_wp2[0]);

  ptsy.push_back(next_wp0[1]);
  ptsy.push_back(next_wp1[1]);
  ptsy.push_back(next_wp2[1]);

  // totally we have five points in ptsx, ptsy

  for(int i = 0; i < ptsx.size(); i++){
    // shift car reference angle to 0 degrees
    double shift_x = ptsx[i] - ref_x;
    double shift_y = ptsy[i] - ref_y;

    ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
    ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
  }

  //create a spline
  tk::spline s;

  // set (x,y) points to the spline
  s.set_points(ptsx, ptsy);

  vector<double> next_x_vals;
  vector<double> next_y_vals;

  // start with all of the previous path points from last time
  for(int i = 0; i < previous_path_y.size(); i++){
    next_x_vals.push_back(previous_path_x[i]);
    next_y_vals.push_back(previous_path_y[i]);
  }

  // calculate how to break up spline points so that we travel at our
  // desired reference speed
  double target_x = 50.0;
  double target_y = s(target_x);
  double target_dist = sqrt(target_x*target_x + target_y*target_y);

  double x_add_on = 0; // we start at local coordinate 0.

  // fill up the rest of our path planner after filling it with previous
  // points
  // double N = (target_dist / 0.02 / ref_vel / 2.237); // 2.237 to convert mph to m/s
  // double delta_x = target_x / N;
  for(int i = 1; i <= 50-previous_path_x.size(); i++){
    double N = (target_dist / (.02*ref_vel) * 2.237);
    double x_point = x_add_on + (target_x)/N;//delta_x;
    double y_point = s(x_point);

    x_add_on = x_point;

    double x_ref = x_point;
    double y_ref = y_point;

    // rotate back to world coordinate
    x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
    y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

    x_point += ref_x;
    y_point += ref_y;

    next_x_vals.push_back(x_point);
    next_y_vals.push_back(y_point);

  }

  return {next_x_vals, next_y_vals};
}


double compute_cost(vector<double> next_x_vals, vector<double> next_y_vals){
  // compute cost of speed with speed limit
  double cost = 0;

  for(int i=1; i < next_x_vals.size(); i++){
    double vx = (next_x_vals[i] - next_x_vals[i-1]) / 0.02;
    double vy = (next_y_vals[i] - next_y_vals[i-1]) / 0.02;

    double cost_i = 49.5 / 2.24 - sqrt(vx*vx + vy*vy);
    cost += cost_i;

  }
  // normalize it over all speed 0
  cost = cost / (next_x_vals.size() - 1) / (49.5 / 2.24);
  // cout << cost << endl;
  return cost;
}


int get_current_lane(double d){
  //lane_number calculation
  // from center line to the most right, lane 1, 2, 3.
  int lane_num;

  if(d < 4 && d >= 0){
    lane_num = 0;
  }
  else if(d < 8 && d >= 4){
    lane_num = 1;
  }
  else if(d < 12 && d >= 8){
    lane_num = 2;
  }
  return lane_num;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  // start on lane 1, the middle lane.
  int lane = 1;
  int target_lane = lane;


  h.onMessage([&lane, &target_lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;


    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;



            int prev_size = previous_path_x.size();

            // sensor Fusion
            // data in sensor_fusion [car's unique ID, car's x position in map coordinates,
            // car's y position in map coordinates, car's x velocity in m/s,
            // car's y velocity in m/s, car's s position in frenet coordinates,
            // car's d position in frenet coordinates] i.e. [id, x, y, vx, vy, s, d]
            if(prev_size > 0){
              car_s = end_path_s;
            }

            bool too_close = false;
            lane = get_current_lane(car_d);


            // below are computed all for the end of previous path, not current time
            // find out the nearest front car index and rear car index of ego vehicle on adjacent lanes
            // find out the nearest front car on ego vehicle's lane

            int front_nearest_car= -1;
            int left_front_nearest_car= -1;
            int right_front_nearest_car= -1;
            int left_rear_nearest_car= -1;
            int right_rear_nearest_car= -1;

            float front_min_dist = 10000;
            float left_front_min_dist = 10000;
            float left_rear_min_dist = 10000;
            float right_front_min_dist = 10000;
            float right_rear_min_dist = 10000;

            for(int i = 0; i < sensor_fusion.size(); i ++){
              vector<double> car_i_data = sensor_fusion[i];
              float d = car_i_data[6];
              double vx = car_i_data[3];
              double vy = car_i_data[4];
              double check_speed = sqrt(vx*vx + vy*vy);
              double check_car_s = car_i_data[5];

              // if using previous value can project s value out
              check_car_s += ((double)prev_size*.02*check_speed);

              double dist = check_car_s - car_s;
              // check my lane
              if (d < (4*(lane+1)) && d > (4*(lane))){

                if((dist > 0) && (dist < front_min_dist)){
                  front_min_dist = dist;
                  front_nearest_car = i;
                }
              }
              // check my left lane
              else if (d < (4*(lane)) && d > (4*(lane-1))){

                // compute left front nearest car
                if((dist > 0) && (dist < left_front_min_dist)){
                  left_front_min_dist = dist;
                  left_front_nearest_car = i;
                }

                // compute left rear nearest car
                if((dist < 0) && (-dist < left_rear_min_dist)){
                  left_rear_min_dist = -dist;
                  left_rear_nearest_car = i;
                }

              }
              // check my right lane
              else if (d < (4*(lane+2)) && d > (4*(lane+1))){

                // compute right front nearest car
                if((dist > 0) && (dist < right_front_min_dist)){
                  right_front_min_dist = dist;
                  right_front_nearest_car = i;
                }

                // compute left rear nearest car
                if((dist < 0) && (-dist < right_rear_min_dist)){
                  right_rear_min_dist = -dist;
                  right_rear_nearest_car = i;
                }

              }
            }
            cout << " " << endl;
            cout << "left: front "<<left_front_nearest_car << " " <<left_front_min_dist <<endl;
            cout << "left: rear " << left_rear_nearest_car << " " << left_rear_min_dist << endl;
            cout << "right: front "<<right_front_nearest_car << " " << right_front_min_dist << endl;
            cout << "right: rear "<< right_rear_nearest_car << " " << right_rear_min_dist << endl;
            // now we know cars around the ego car
            // if left_front_nearest_car[0] == -1, means left_front_nearest_car not exist

            vector<vector<double>> next_vals;
            vector<double> next_x_vals;
            vector<double> next_y_vals;

            double cost_lk;
            double cost_lcl = 10000;
            double cost_lcr = 10000;


            // behavior planner
            // we perform behavior planning if target_lane equals lane.
            // if not equal, means we are changing lanes.
            if(target_lane == lane){

              if (front_nearest_car != -1){
                // if there are front vehicles

                // check whether need to keep lane
                if (front_min_dist < 30){
                  // if front vehicle is less than 30m to ego vehicle
                  too_close = true;

                  // decide keep lane, or lane change left, or lane change right

                  // lane keep cost
                  next_vals  = generate_path(too_close, lane, j[1],
                    map_waypoints_x, map_waypoints_y, map_waypoints_s);

                  next_x_vals = next_vals[0];
                  next_y_vals = next_vals[1];


                  cost_lk = compute_cost(next_x_vals, next_y_vals);

                  too_close = true;


                  bool able_to_change;
                  // lanes 1 and 2, compute change left cost
                  if(lane > 0){
                    able_to_change = false;

                    if(left_front_nearest_car == -1 && (left_rear_nearest_car == -1 ||
                    left_rear_min_dist > 10)){
                      // no front car, able to change
                      able_to_change = true;

                    }else if((left_front_nearest_car != -1 && left_front_min_dist > 40) && (left_rear_nearest_car == -1 ||
                      left_rear_min_dist > 10)){
                        // if front nearest car is more than 50m away from me and rear car is 10m
                        // away from me
                        able_to_change = true;

                    }

                    if(able_to_change){
                      next_vals  = generate_path(too_close, lane-1, j[1],
                        map_waypoints_x, map_waypoints_y, map_waypoints_s);

                      next_x_vals = next_vals[0];
                      next_y_vals = next_vals[1];

                      cost_lcl = compute_cost(next_x_vals, next_y_vals);
                    }
                  // end of lane change left cost computation
                  }


                  // lanes 0 and 1, compute change right cost
                  if(lane < 2){
                    // defaultly, give large cost for not able to change lane.
                    able_to_change = false;

                    if(right_front_nearest_car == -1 && (right_rear_nearest_car == -1 || right_rear_min_dist > 10)){
                      // if there is no front car on right lane
                      able_to_change = true;

                    }else if((right_front_nearest_car != -1 && right_front_min_dist > 40)
                              && (right_rear_nearest_car == -1 || right_rear_min_dist > 10)){
                        // if rear car exist, it should leave some space for ego car
                        able_to_change = true;
                    }


                    if(able_to_change){
                      // if it is able to change lane right, compute it.
                      next_vals  = generate_path(too_close, lane+1, j[1],
                    map_waypoints_x, map_waypoints_y, map_waypoints_s);

                      next_x_vals = next_vals[0];
                      next_y_vals = next_vals[1];

                      cost_lcr = compute_cost(next_x_vals, next_y_vals);
                    }
                  // end of lane change right cost computation
                  }


                  /////// choose the trajectory of minimum cost
                  // double costs [] = {cost_lk, cost_lcl, cost_lcr};
                  // should have some simple functions to use, but I didn't find
                  double min_cost = 1000;
                  int min_cost_idx;

                  if(cost_lcl == 10000 && cost_lcr == 10000){
                    // both change left and right are not possible
                    target_lane = lane;
                  }else{
                    double costs [] = {cost_lcl, cost_lcr};
                    min_cost = std::min(cost_lcl, cost_lcr);
                    // cout<< "costs: "<<cost_lcl << "  "<<cost_lcr<< "  cost min" << min_cost <<endl;

                    if(min_cost == cost_lcl){
                      target_lane = lane - 1;
                    }else{
                      target_lane = lane + 1;
                    }
                  }

                  // if(min_cost_idx == 1){
                  //   target_lane = lane - 1;
                  // }else if(min_cost_idx==2){
                  //   target_lane = lane + 1;
                  // }else{
                  //   target_lane = lane;
                  // }
                  cout<<"               minimum cost is " << min_cost <<endl;
                  cout<<"               selected lane " << target_lane <<endl;


                    ////// end of choosing trajectory of minimum cost

                }
                else{
                  // front nearest vehicle is far away, keep lane
                  too_close = false;
                  target_lane = lane;
                }
              }else{
                // no front vehicle, keep lane
                too_close = false;
                target_lane = lane;
              }


            }else{
              // changing lane
              too_close = true;
            }

            // correct weird behavior of left changing lane on lane 0
            if(target_lane < 0){
              target_lane = 0;
            }else if(target_lane > 2){
              target_lane = 2;
            }

            next_vals  = generate_path(too_close, target_lane, j[1],
              map_waypoints_x, map_waypoints_y, map_waypoints_s);
            next_x_vals = next_vals[0];
            next_y_vals = next_vals[1];


          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
