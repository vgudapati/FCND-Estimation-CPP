#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;

#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();

  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);

  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to
  //   individual motor thrust commands
  // INPUTS:
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS:
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  /*
  cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right
   */

  /*
  Equations for generating the desiredThrusts

  We know;

  tau_x = (F1 - F2 - F3 + F4) * l
  tau_y = (F1 + F2 - F3 - F4) * l
  tau_z = tau_1 + tau_2 + tau_3 + tau_4
  F = F1 + F2 + F3 + F4

   We need to solve these for F1, F2, F3, F4 (desired thrusts)

  We also know;

  tau_1 = k_m * omega1 ** 2;     F1 = k_f * omega1 ** 2  =>  tau_1 = k_m * (F1 / k_f)
  tau_2 = - k_m * omega2 ** 2;   F2 = k_f * omega2 ** 2  =>  tau_2 = -k_m * (F2 / k_f)
  tau_3 = k_m * omega3 ** 2;     F3 = k_f * omega3 ** 2  =>  tau_3 = k_m * (F3 / k_f)
  tau_4 = -k_m * omega ** 4;     F4 = k_f * omega4 ** 2  =>  tau_4 = -k_m * (F4 / k_f)

  This leads to;

  tau_z = kappa * (F1 - F2 + F3 - F4), where kappa is k_m/k_f (drag/thrust)
  l = L / (2 * sqrt(2))

  This gives us

  F1 - F2 - F3 + F4 = tau_x/l                       -       (1)
  F1 + F2 - F3 - F4 = tau_y/l                       -       (2)
  F1 - F2 + F3 - F4 = tau_z / kappa                 -       (3)
  F1 + F2 + F3 + F4 = F                             -       (4)

   (1) + (4) gives us;

   2 * F1 + 2 * F4 = (tau_x/l) + F                  -       (5)

   (2) + (3) gives us;

   2 * F1 - 2 * F4 = (tau_y/l) + (tau_z / kappa)    -       (6)

   (4) - (1) gives us;

   2 * F2 + 2 * F3 = F - tau_x / l;                 -       (7)

   (2) - (3) gives us;

   2 * F2 - 2 * F3 = tau_y/l - tau_z / kappa        -       (8)

   (5) + (6) gives us;

   4 * F1 = (tau_x/l) + F + (tau_y/l) + (tau_z / kappa) =>  F1 = F/4 + (tau_x)/l/4 + tau_y/l/4 + tau_z/kappa/4

   (7) + (8) gives us;

   4 * F2 = F - tau_x / l + tau_y / l - tau_z / kappa   =>  F2 = F/4 - (tau_x)/l/4 + tau_y/l/4 - tau_z/kappa/4

   (7) - (8) gives us;

   4 * F3 = F - tau_x / l - tau_y/l + tau_z / kappa     =>  F3 = F/4 - tau_x/l/4 - tau_y/l/4 + tau_z/kappa/4

   (5) - (6) gives us;

   4 * F4 = F + tau_x/l - tau_y/l - tau_z / kappa       =>  F4 = F/4 + tau_x/l/4 - tau_y/l/4 - tau_z/kappa/4

   Since tau's are moments and F is the total force, let us relate them to inputs

   tau_x = momentCmd.x
   tau_y = momentCmd.y
   tau_z = -momentCmd.z
   F = collThrustCmd

   Let us code with these values and the above equations and translating l to L.

  */
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  float F1 = collThrustCmd/4.f + momentCmd.x/L /2.f/ sqrtf(2) + momentCmd.y/L /2.f/ sqrtf(2) - momentCmd.z/kappa/4.f;
  float F2 = collThrustCmd/4.f - momentCmd.x/L /2.f/ sqrtf(2) + momentCmd.y/L /2.f/ sqrtf(2) + momentCmd.z/kappa/4.f;
  float F3 = collThrustCmd/4.f - momentCmd.x/L /2.f/ sqrtf(2) - momentCmd.y/L /2.f/ sqrtf(2) - momentCmd.z/kappa/4.f;
  float F4 = collThrustCmd/4.f + momentCmd.x/L /2.f/ sqrtf(2) - momentCmd.y/L /2.f/ sqrtf(2) + momentCmd.z/kappa/4.f;

  cmd.desiredThrustsN[0] = F1;    //  front left
  cmd.desiredThrustsN[1] = F2;    //  front right
  cmd.desiredThrustsN[2] = F4;    //  rear left
  cmd.desiredThrustsN[3] = F3;    //  rear right

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS:
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS:
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
  We wll be coding using equations

  p_error = p_target - p_actual
  u_bar_p = k_p_p * p_error
  and tau_x = Ix * u_bar_p

  */

  V3F p_error = pqrCmd - pqr;
  V3F omega_dot_desired = kpPQR * p_error;
  momentCmd = V3F(Ixx, Iyy, Izz) * omega_dot_desired;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS:
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  /*

   We have the formulas;

   b_x_c_dot = k_p * ( b_x_c - b_x_a), where b_x_a = R13
   b_y_c_dot = k_p * ( b_y_c - b_y_a), where b_y_a = R23

   Also, the angular velocities (p_c, q_c) = (1/R33) * ((R21, -R11), (R22, -R12)) * (b_x_c_dot, b_y_c_dot)

   This gives us;

   p_c = (1/R33) * (R21 * b_x_c_dot - R11 * b_y_c_dot)

   =>

   p_c = (1/R33) * (R21 * k_p * (b_x_c - R13) - R11 * k_p * (b_y_c - R23))

   q_c = (1/R33) * (R22 * b_x_c_dot - R12 * b_y_c_dot)

   =>

   q_c = (1/R33) * (R22 * k_p * (b_x_c - R13) - R12 * k_p * (b_y_c - R23))

   We also know;

   c = F/m

   p_c and q_c are the roll and putch rates

   So, pqrCmd.x = p_c; pqrCmd.y = q_c; and we will set pqrCmd.z = 0;

   With this we can now code
  */

  float b_x_c = -CONSTRAIN(accelCmd[0] / (collThrustCmd / mass), -maxTiltAngle, maxTiltAngle);
  float b_y_c = -CONSTRAIN(accelCmd[1] / (collThrustCmd / mass), -maxTiltAngle, maxTiltAngle);

  if (collThrustCmd < 0 ) {
    b_x_c = 0;
    b_y_c = 0;
  }

  pqrCmd.x = (1 / R(2, 2))*(-R(1, 0) * kpBank*(R(0, 2) - b_x_c) + R(0, 0) * kpBank*(R(1, 2) - b_y_c));
  pqrCmd.y = (1 / R(2, 2))*(-R(1, 1) * kpBank*(R(0, 2) - b_x_c) + R(0, 1) * kpBank*(R(1, 2) - b_y_c));
  pqrCmd.z = 0;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical
  //   acceleration feed-forward command
  // INPUTS:
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*

    We have formulas;

    c = (u_1_bar - g)/b_z, where b_z = R33;
    u_1_bar = k_p_z * (z_target - z_atual) * k_d_z(z_dot_target - z_dot_zctual) + z_dot_dot_target
    and F = c * m

    using these formulas and integrating the altitude error we can code as below
  */

  velZCmd += kpPosZ * (posZCmd - posZ);
  integratedAltitudeError += (posZCmd - posZ) * dt;
  velZCmd = CONSTRAIN(velZCmd, -maxAscentRate, maxDescentRate);

  thrust = - ((kpVelZ * (velZCmd - velZ) + KiPosZ * integratedAltitudeError + accelZCmd - 9.81f)/R(2,2) * mass);


  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS:
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations.
  //     the Z component should be 0
  // HINTS:
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////


  /*

   From the python implementation in the excercises we have

   x_dot_dot_command = x_k_p * (x_target - x_actual) + x_k_d * (x_dot_target - x_dot_actual) + x_dot_dot_target

   realigning the inputs here as:

   x_dot_dot_target = accelCmdFF;
   x_target = posCmd
   x_actual = pos
   x_dot_target = velCmd
   x_dot_actual = vel

   and limiting the maximum horizontal velocity and acceleration to maxSpeedXY and maxAccelXY as per requirements, we can code as below.

  */
  velCmd += kpPosXY * (posCmd - pos);

  if (velCmd.mag() > maxSpeedXY)
  {
    velCmd = velCmd * maxSpeedXY / velCmd.mag();
  }

  accelCmd += kpVelXY * (velCmd - vel);
  if (accelCmd.mag() > maxAccelXY)
  {
    accelCmd = accelCmd * maxAccelXY / accelCmd.mag();
  }


  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS:
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS:
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b].
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
    Using the equation;

    r_c = k_p_yaw * ( psi_target - psi_actual),

    we can code as below.

  */

  float yawError = yawCmd - yaw;
  yawError = fmodf(yawError, F_PI*2.f);
  if (yawError > F_PI)
  {
    yawError -= 2.f * F_PI;
  }
  else if (yawError < -F_PI)
  {
    yawError += 2.f * F_PI;
  }
  yawRateCmd = yawError * kpYaw;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);

  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);

  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
