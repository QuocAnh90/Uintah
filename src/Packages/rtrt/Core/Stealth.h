

#ifndef STEALTH_H
#define STEALTH_H 1

#include <iostream>
#include <vector>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

namespace rtrt {

using namespace SCIRun;

using std::cout;
using std::vector;

// "Stealth" comes from the DoD term for an invisible watcher on the
// simulated battlefield.  The stealths are used in relation to
// cameras.  They contain movement information, but the camera
// contains direction information.
  
class Stealth {

  ////////////////////////////////////////////////////////
  // Stealth Movement
  //
  //   This camera does NOT move like an airplane.
  // The driver (user) can specify movement in the following
  // manner.
  //
  //     +  Accelerate Forward
  //     -  Accelerate Backward
  //     <- Begin turning left 
  //     -> Begin turning right
  //     ^  Begin pitching forward (look down)
  //     v  Begin pitching backward (look up)
  //     7 (keypad) Accelerate to the left
  //     9 (keypad) Accelerate to the right
  //     0 (keypad) Stop ALL movement
  //     . (keypad) Slow down (in all movements (including turns))
  //
  //   Acceleration is in respect to the direction that the eye is
  // looking.  

public:

  // Scale of '50' is good if the eye is 5 units (or so) from
  // the area of interest and the area of interest is ~20 units across.
  Stealth( double scale, double gravity_force );
  ~Stealth();

  // Tells the eye (camera) to update its position based on its current
  // velocity vector.
  void updatePosition();

  inline double getSpeed( int direction ) const {
    switch( direction ) {
    case 0:
      return speed_;
    case 1:
      return horizontal_speed_;
    case 2:
      return vertical_speed_;
    case 3:
      return pitch_speed_;
    case 4:
      return rotate_speed_;
    default:
      cout << "Error in Stealth::getSpeed, bad direction " << direction 
	   << "\n";
      exit( 1 );
      return 0.0;
    }
  }

  // Slows down in all dimensions (pitch, speed, turn, etc);
  void slowDown();
  void stopAllMovement();
  void stopPitchAndRotate();
  void stopPitch();
  void stopRotate();

  void slideLeft();
  void slideRight();

  void goUp();
  void goDown();

  void accelerate();
  void decelerate();
  void turnRight();
  void turnLeft();
  void pitchUp();
  void pitchDown();

  // Display the Stealth's speeds, etc.
  void print();

  // Returns next location in the path and the new view vector.
  void getNextLocation( Point & point, Point & look_at );

  // Clear out path stealth is to follow.
  void clearPath();
  void addToPath( const Point & eye, const Point & look_at );
  void loadPath( const string & filename );
  void savePath( const string & filename );

  // If gravity is on, the stealth/camera will "fall" towards the ground.
  void toggleGravity();

  // Stealths and Cameras are highly integrated right now... perhaps
  // stealth should be a sub class in camera?  This needs more thought.

  bool   gravityIsOn() { return gravity_on_; }
  double getGravityForce() { return gravity_force_; }

  bool   moving(); // Returns true if moving in any direction

private:

  void increase_a_speed( double & speed, int & accel_cnt );
  void decrease_a_speed( double & speed, int & accel_cnt );

  // Scale is based on the size of the "universe".  It effects how fast
  // the stealth will move.
  double scale_;

  // Speeds (in units per frame)
  double speed_;            // + forward, - backward
  double horizontal_speed_; // + right, - left
  double vertical_speed_;   // + up, - down
  double pitch_speed_;      // + down, - up
  double rotate_speed_;     // + right, - left

  // Acceleration counts represent the number of velocity change 
  // requests that the user has made.  They are used to judge how
  // much velocity to add to the velocity vector.  The more requests,
  // the faster we "accelerate" for the next similar request.
  int accel_cnt_;            // + forward, - backward. 
  int horizontal_accel_cnt_; // + right, - left
  int vertical_accel_cnt_;   // + up, - down
  int pitch_accel_cnt_;      // + up, - down
  int rotate_accel_cnt_;     // + right, - left

  // Path information  (rough draft)

  vector<Point> path_;     // Points along the path for the stealth to move.
  vector<Point> look_ats_; // Look at locations along path.

  double        segment_percentage_;

  bool          gravity_on_;
  double        gravity_force_;
};

} // end namespace rtrt

#endif
