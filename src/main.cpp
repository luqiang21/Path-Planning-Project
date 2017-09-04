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

int getLaneNumber(double d){
  // lane_number calculation
  // from center line to the most right, lane 1, 2, 3.
  int lane_num;

  if(d < 4 && d > 0){
    lane_num = 1;
  }
  else if(d < 8 && d > 4){
    lane_num = 2;
  }else if(d < 12 && d > 8){
    lane_num = 3;
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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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



          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

             double pos_x;
             double pos_y;
             double angle;
             int path_size = previous_path_x.size();

            //  for(int i = 0; i < path_size; i++)
            //  {
            //      next_x_vals.push_back(previous_path_x[i]);
            //      next_y_vals.push_back(previous_path_y[i]);
            //  }


            if(path_size != 0)
          {
              next_x_vals.push_back(previous_path_x[path_size-2]);
              next_y_vals.push_back(previous_path_y[path_size-2]);
              next_x_vals.push_back(previous_path_x[path_size-1]);
              next_y_vals.push_back(previous_path_y[path_size-1]);
            }

              if(path_size == 0)
            {
                pos_x = car_x;
                pos_y = car_y;
                angle = deg2rad(car_yaw);
            }
            else
            {
                pos_x = previous_path_x[path_size-1];
                pos_y = previous_path_y[path_size-1];

                double pos_x2 = previous_path_x[path_size-2];
                double pos_y2 = previous_path_y[path_size-2];
                angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
            }

             // // // keep on its current lane
             // find out the front vehicle of the ego vehicle
             cout << sensor_fusion << endl;
            int car_lane_num = getLaneNumber(car_d);
            cout << " my d is "<< car_d << "lane is " << car_lane_num << endl;

            int lane_num = -1;
            int front_car_id = -1;
            double front_car_dist, car_length = 5.; // need to verify the length.
            int front_car_lane = -1;
            double front_car_vx, front_car_vy;

            for(int i = 0; i < sensor_fusion.size(); i++)
            {
                auto sensor_data  = sensor_fusion[i];
                int id = sensor_data[0];
                double x = sensor_data[1], y = sensor_data[2];
                double vx = sensor_data[3], vy = sensor_data[4];
                double s = sensor_data[5], d = sensor_data[6];

                lane_num = getLaneNumber(d);

                cout << id << "   " << lane_num << endl;
                cout << "d " << d << endl;
                if(lane_num == car_lane_num){
                  double dist = distance(car_x, car_y, x, y) - car_length;

                  if(front_car_id == -1 || dist < front_car_dist){
                    front_car_id = id;
                    front_car_dist = dist;
                    front_car_lane = lane_num;
                    front_car_vx = vx;
                    front_car_vy = vy;
                  }
                }
            }
            cout << "the immediate front car of me is " << front_car_id << " lane " << front_car_lane<<" distance " << front_car_dist<< endl;

             // obtain the closest way point.
            //  pos_x = previous_path_x[path_size-1];
            //  pos_y = previous_path_y[path_size-1];

            //  int closet_way_point_idx = ClosestWaypoint(pos_x, pos_y, map_waypoints_x, map_waypoints_y);
             //
            //  double s = map_waypoints_s[closet_way_point_idx];
            //  double d_x = map_waypoints_dx[closet_way_point_idx];
            //  double d_y = map_waypoints_dy[closet_way_point_idx];
             //
            //  double d_magnitude = sqrt(d_x * d_x + d_y * d_y);
            double s = car_s;
            //  if( (d_magnitude - 2) <=
            double dist_inc = 0.4; // distance increased in 1/50 seconds.

            // set ego car's speed based on front car's speed using car following model.

            double front_car_v = sqrt(front_car_vx*front_car_vx + front_car_vy*front_car_vy);
            front_car_v = front_car_v / 2.237; // convert mph to m/s

            cout << "car speed is "<<car_speed <<endl;
            cout << "immediate front car's speed is " << front_car_v << endl;
            double safe_v = sqrt(front_car_v*front_car_v + 2*5*(front_car_dist - 5)); // safe minimum gap 2m
            if(front_car_id == -1){
              safe_v = 25.;
            }

            car_speed = car_speed / 2.237 ; // convert mph to m/s
            if(car_speed > safe_v){
              // if the ego car exceeds safe speed, set to safe speed.
              dist_inc = safe_v / 50;
              cout << "vehicle's speed is set to safe speed " << safe_v << endl;
              cout << "dist_inc is " << dist_inc << endl;
            }



            double max_decel = 5.;
            for(int i = 0; i < 50 - 2; i++)
            {
                if(car_speed > safe_v){
                  car_speed -= max_decel * 1/50.;

                }
                if(car_speed > 25.){
                  car_speed = 25;

                }
                if(car_speed < safe_v){
                  car_speed = safe_v;
                }
                cout << "car speed " << car_speed << endl;
                dist_inc = car_speed / 50.;
                cout << "safe speed " << safe_v <<endl;
                cout << "dist_inc "<< dist_inc <<endl;

                s += dist_inc;

                auto x_y = getXY(s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

                // cout << x_y << endl;
                double x = x_y[0];
                double y = x_y[1];
                next_x_vals.push_back(x);
                next_y_vals.push_back(y);

            }
            cout << "hahaha"<< endl;



            //  NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
             // compute frenet coordinates by assuming a constant speed

             // convert frenet coordinates to x, y


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
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
