
#include <Packages/rtrt/Core/Stealth.h>

#include <math.h>

using namespace rtrt;
using namespace SCIRun;
using std::cerr;

Stealth::Stealth( double scale /* = 100 */ ) :
  speed_( 0 ), horizontal_speed_( 0 ), vertical_speed_( 0 ),
  pitch_speed_( 0 ), rotate_speed_( 0 ),
  accel_cnt_( 0 ), horizontal_accel_cnt_( 0 ),
  vertical_accel_cnt_( 0 ), pitch_accel_cnt_( 0 ), rotate_accel_cnt_( 0 ),
  scale_( scale ), segment_percentage_( 0 )
{
  path_.push_back( Point(-5,-5,5) );
  path_.push_back( Point(5,5,1) );
  path_.push_back( Point(0,10,5) );
  path_.push_back( Point(5,5,1) );
  path_.push_back( Point(-5,-5,5) );
}

Stealth::~Stealth()
{
}

void
Stealth::print()
{
  cerr << "Stealth print not done yet!\n";
}

void
Stealth::slideRight()
{
  increase_a_speed( horizontal_speed_, horizontal_accel_cnt_ );
}

void
Stealth::slideLeft()
{
  decrease_a_speed( horizontal_speed_, horizontal_accel_cnt_ );
}

void
Stealth::goUp()
{
  increase_a_speed( vertical_speed_, vertical_accel_cnt_ );
}

void
Stealth::goDown()
{
  decrease_a_speed( vertical_speed_, vertical_accel_cnt_ );
}

void
Stealth::accelerate()
{
  increase_a_speed( speed_, accel_cnt_ );
}

void
Stealth::decelerate()
{
  decrease_a_speed( speed_, accel_cnt_ );
}

void
Stealth::increase_a_speed( double & speed, int & accel_cnt )
{
  accel_cnt++;

  // Amount to scale acceleration by.  Probably need to find out "size"
  // (dimensions?) of data set to set this appropriately.
  double scale = 100;
  
  // Amount of acceleration doubles for each consecutive request.
  double amount;

  if( accel_cnt == 0 )
    amount = 1 / scale;
  else if( accel_cnt > 0 )
    amount = pow(2.0, (double)accel_cnt - 1.0) / scale;
  else
    amount = pow(2.0, (double)(-accel_cnt)) / scale;

  printf("accelerating by %lf\n", amount);

  speed += amount;
}

void
Stealth::decrease_a_speed( double & speed, int & accel_cnt )
{
  accel_cnt--;

  // Amount to scale acceleration by.  Probably need to find out "size"
  // (dimensions?) of data set to set this appropriately.
  double scale = 100;
  
  // Amount of acceleration doubles for each consecutive request.
  double amount = pow(2.0, (double)abs(accel_cnt)) / scale;

  if( accel_cnt >= 0 )
    amount = pow(2.0, accel_cnt) / scale;
  else
    amount = pow(2.0, (double)(-accel_cnt) - 1.0) / scale;

  printf("reversing by %lf\n", amount);

  speed -= amount;
}

void
Stealth::turnRight()
{
  increase_a_speed( rotate_speed_, rotate_accel_cnt_ );
}

void
Stealth::turnLeft()
{
  decrease_a_speed( rotate_speed_, rotate_accel_cnt_ );
}

void
Stealth::pitchUp()
{
  increase_a_speed( pitch_speed_, pitch_accel_cnt_ );
}

void
Stealth::pitchDown()
{
  decrease_a_speed( pitch_speed_, pitch_accel_cnt_ );
}

void
Stealth::stopAllMovement()
{
  speed_ = 0;
  horizontal_speed_ = 0;
  vertical_speed_ = 0;
  pitch_speed_ = 0;
  rotate_speed_ = 0;

  accel_cnt_ = 0;
  horizontal_accel_cnt_ = 0;
  vertical_accel_cnt_ = 0;
  pitch_accel_cnt_ = 0;
  rotate_accel_cnt_ = 0;
}

void
Stealth::slowDown()
{
  double slow_factor = .9;
  speed_            /= slow_factor;
  horizontal_speed_ /= slow_factor;
  vertical_speed_   /= slow_factor;
  pitch_speed_      /= slow_factor;
  rotate_speed_     /= slow_factor;
}

Point
Stealth::getNextLocation()
{
  if( path_.size() < 2 ) {
    cout << "Not enough path enformation!\n";
    return Point(0,0,0);
  }

  segment_percentage_ += 1;

  int begin_index = segment_percentage_ / 100;
  int end_index   = (segment_percentage_ / 100) + 1;

  if( end_index == path_.size() )
    {
      cout << "got to the end of the path... restarting\n";
      segment_percentage_ = 1;
      begin_index = 0; end_index = 1;
    }

  Point & begin = path_[ begin_index ];
  Point & end = path_[ end_index ];

  double local_percent = segment_percentage_ - (100*begin_index);

  Point new_loc = begin + ((end - begin) * (local_percent / 100.0));

  return new_loc;
}
