/*
 Main function done as part of diploma thesis
 Data of last edition: 10.06.2019

 */

#include <unistd.h>
#include <MPU9250_Master_I2C.h>  // do MPU9250 po i2c
#include <stdint.h>
#include <stdio.h>
#include <wiringPi.h>
#include <wiringSerial.h>
#include <errmsg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include <eigen/Eigen/Dense>

#include <pthread.h>

#include <pigpio.h>
#include <stdlib.h>
#include <mcp3004.h>           //do MCP3008      

#include <eigen/Eigen/Dense>

#include "uNavAHRS.h" 
#include "mahonyMadgwcickFilter.h"

extern "C" {
#include  "kalman.h"
}Kalm anFilter f = alloc_filter(4, 4);

pthread_t watekImuEuler;
pthread_t watekImuReadData;  // odczyt z maksymalna czestoltiwoscia
pthread_t watekGpsReadData;
pthread_t watekADCReadData;
pthread_t watekRadioData;
pthread_t watekMotorsWrite;
pthread_t watekPostion;
pthread_t watekLogData;
pthread_t watekEKF;

//do MPU
float ax, ay, az, gx, gy, gz, mx, my, mz, temperature;

float axb, ayb, azb, gxb, gyb, gzb;

// M M FILTERS
float q[4];
float beta = 0.1;
float deltat = 1e-3;
float Kp = 1;
float Ki = 1;
float eInt[4];

Eigen::Matrix<float, 4, 1> q2;
Eigen::Matrix<float, 3, 1> eul;

//float bx,by,bz,sx,sy,sz;

//premeasured
float bx = 82.442993;
float by = 333.731293;
float bz = -88.417122;
float sx = 0.964680;
float sy = 1.087065;
float sz = 0.958333;

static const MPUIMU::Gscale_t GSCALE = MPUIMU::GFS_250DPS;
static const MPUIMU::Ascale_t ASCALE = MPUIMU::AFS_2G;
static const MPU9250::Mscale_t MSCALE = MPU9250::MFS_16BITS;
static const MPU9250::Mmode_t MMODE = MPU9250::M_100Hz;
static const uint8_t SAMPLE_RATE_DIVISOR = 0x02;

//od MCP
#define BASE 200
#define SPI_CHAN 0

//gpio zgdone z BCM do sterowania silnikami DC

#define MOTOR_A_FORWARD 12
#define MOTOR_A_BACKWARD 16

#define MOTOR_B_FORWARD 27
#define MOTOR_B_BACKWARD 22

//offsety czujnikow pradu

#define prad_motor_A_offset 80
#define prad_motor_B_offset 190

using Eigen::MatrixXd;

//acc biases
double acc_bias_x;  //acc biases
double acc_bias_y;
double acc_bias_z;

//gyro biases
double gyro_bias_x;
double gyro_bias_y;
double gyro_bias_z;

//gyro yaw correction
double gyroYawFactor = 1;

//do poprawy  na coĹ› bardziej zanaczacego
int motorsPrint;  // motors
int eulerPrint;
int velAndPosPrint;  //
int senosorsRawPrint = 0;
int printEKF = 0;

///////////// dane bez filtracji
float pitchAcc;
float rollAcc;
float yawMag;

float pitchGyro;
float yawGyro;
float rollGyro;

int yaw_dest = 45;
int aRead[8]; // odczyty z MCP3008

bool offsetCaluclatedAcc = false;
bool offsetCaluclatedGyro = false;
//prady
float prad_motor_A = 0;
float suma_prad_motor_A = 0;

float prad_motor_B = 0;
float suma_prad_motor_B = 0;

//radio data variables
char data[50];
int daneRx[8];
int fd_Radio;

// filtr EKF

uNavAHRS filtrEKF;

//do switch case
int polecenie = 0;

//Struktura glowna
struct Boat {
	int motorsArmed;

	int speed_motor_A;
	int speed_motor_B;

	float yaw;
	float pitch;
	float roll;

	float x_pos;
	float y_pos;

	float x_base;
	float y_base;

} suv;

struct anglesEKF {
	float yaw;
	float pitch;
	float roll;
} ekf;

struct RadioChanels {  //kolejnosc odbierania
	//0 -ramka 1155
	uint16_t armed;           //1
	uint16_t left_right;      //2
	uint16_t forward_backward;      //3
	uint16_t set_point;       //4
	uint16_t send_gps_pos;    //5
	uint16_t auto_mode;       //6

} RadioChannel;

// Pin definitions for MPU9250
static const uint8_t intPin = 0;   //  MPU9250 interrupt

// Interrupt support
static bool gotNewData;

static void myinthandler() {
	gotNewData = true;
}

// Instantiate MPU9250 class in master mode
static MPU9250_Master_I2C imu(ASCALE, GSCALE, MSCALE, MMODE,
		SAMPLE_RATE_DIVISOR);

int fd;

Eigen::Matrix<float, 4, 1> Eul2Quat(Eigen::Matrix<float, 3, 1> eul) {
	Eigen::Matrix<float, 4, 1> q;
	q(0, 0) = cos(eul(2, 0) / 2.0) * cos(eul(1, 0) / 2.0) * cos(eul(0, 0) / 2.0)
			+ sin(eul(2, 0) / 2.0) * sin(eul(1, 0) / 2.0)
					* sin(eul(0, 0) / 2.0);
	q(1, 0) = cos(eul(2, 0) / 2.0) * cos(eul(1, 0) / 2.0) * sin(eul(0, 0) / 2.0)
			- sin(eul(2, 0) / 2.0) * sin(eul(1, 0) / 2.0)
					* cos(eul(0, 0) / 2.0);
	q(2, 0) = cos(eul(2, 0) / 2.0) * sin(eul(1, 0) / 2.0) * cos(eul(0, 0) / 2.0)
			+ sin(eul(2, 0) / 2.0) * cos(eul(1, 0) / 2.0)
					* sin(eul(0, 0) / 2.0);
	q(3, 0) = sin(eul(2, 0) / 2.0) * cos(eul(1, 0) / 2.0) * cos(eul(0, 0) / 2.0)
			- cos(eul(2, 0) / 2.0) * sin(eul(1, 0) / 2.0)
					* sin(eul(0, 0) / 2.0);
	return q;
}

Eigen::Matrix<float, 3, 1> Quat2Eul(Eigen::Matrix<float, 4, 1> q) {
	Eigen::Matrix<float, 3, 1> eul;
	eul(0, 0) = atan2(2.0 * q(2, 0) * q(3, 0) + 2.0 * q(0, 0) * q(1, 0),
			2.0 * q(0, 0) * q(0, 0) + 2.0l * q(3, 0) * q(3, 0) - 1.0);
	eul(1, 0) = asin(-2.0 * q(1, 0) * q(3, 0) + 2.0 * q(0, 0) * q(2, 0));
	eul(2, 0) = atan2(2.0 * q(1, 0) * q(2, 0) + 2.0 * q(0, 0) * q(3, 0),
			2.0 * q(0, 0) * q(0, 0) + 2.0l * q(1, 0) * q(1, 0) - 1.0);
	return eul;
}

float reg_PI(float error, float Kp, float Ki, float up_limit, float down_limit,
		float dT) {  // uchyb,Kp, Ki,nasycenia,dT{
	static float integral = 0;
	float output;

	integral += error * dT;
	if (integral > up_limit)
		integral = up_limit;  //anti wind-up for integral part
	if (integral < down_limit)
		integral = down_limit;

	output = (Kp * error) + (Ki * integral);

	if (output > up_limit)
		output = up_limit;  //anti wind-up for final output
	if (output < down_limit)
		output = down_limit;

	return output;

}

float angle360(float dta) {

	if (dta < 0.0)
		dta += 2.0 * M_PI;

	if (dta > (2.0 * M_PI))
		dta -= 2 * 0 * M_PI;

	return dta;
}

void gpsInit() {

	if ((fd = serialOpen("/dev/ttyAMA0", 9600)) < 0) {
		fprintf(stderr, "Unable to open serial device: %s\n", strerror(errno));
		//return 1 ;
	}

	if (wiringPiSetup() == -1) {
		fprintf(stdout, "Unable to start wiringPi: %s\n", strerror(errno));
		// return 1 ;
	}
}

void pigpioInit() {  //inicializacj pigpioInit
	if (gpioInitialise() < 0) {
		fprintf(stderr, "pigpio initialisation failed\n");

	}

}

void *imuRead(void *arg) {
	int delay_time = (int) arg;
	static int x = 0;

	while (1) {
		if (gotNewData) { // On interrupt, read data

			gotNewData = false;     // reset gotNewData flag

			if (imu.checkNewData()) { //odbiĂłr danych z IMU
				x++;
				imu.readAccelerometer(axb, ayb, az);
				imu.readGyrometer(gxb, gyb, gzb);

				if (offsetCaluclatedAcc) {   //uwzglednianie offsetow

					ax = axb - acc_bias_x;
					ay = ayb - acc_bias_y;

					gx = gxb - gyro_bias_x;
					gy = gyb - gyro_bias_y;
					gz = gzb - gyro_bias_z;

					if (senosorsRawPrint && x > 10) {
						x = 0;
						printf("ax = %5d  ay = %5d  az = %5d mg   ",
								(int) (1000 * ax), (int) (1000 * ay),
								(int) (1000 * az));
						printf("gx = %+2.2f  gy = %+2.2f  gz = %+2.2f deg/s\n",
								gx, gy, gz);

					}

				} else {
					ax = axb;
					ay = ayb;
					gx = gxb;
					gy = gyb;
					gz = gzb;
				}

				imu.readMagnetometer(mx, my, mz);
				temperature = imu.readTemperature();

			}

		}

		delay(delay_time);
	}
	return NULL;
}

void calibrateGyroyYaw() {
	float yawInitial = yawGyro;
	int enter;
	printf(
			"Wykonaj peĹ‚ny obrĂłt gyroskpou w prawÄ… stronÄ™ [360 deg] i zatwierdz [1]\n ");
	scanf("%d", &enter);
	if (enter == 1) {
		gyroYawFactor = (-yawGyro + yawInitial) / (2 * M_PI);
		printf(
				"ObrĂłt o 360, wartoĹ›c zmierzona przez gyroskop =  %4.4f, czynik korekcji = %4.4f \n",
				(yawGyro - yawInitial) * 57.3, gyroYawFactor);
		yawGyro = yawInitial;
	}

}

#define NUM_OF_SAMPLES 2048
#define GRAVITY  9.8
// Watek do obliczania pozycji X i Y
void *imuPositionCalc(void *arg) {

	int delay_time;
	delay_time = (int) arg;
	float dt = delay_time / 1000.0;
	int i = 0;
	int j = 0;
	static float velocity[3]; //X, Y, Z

	static float velocity_1[3];
	static float position[3];  //X, Y,Z

	static double ax_avg = 0;
	static double ay_avg = 0;
	static double az_avg = 0;
	//poprzednie probki
	static double ax_1 = 0;
	static double ay_1 = 0;
	//float actual_pos;

	static int count_end_mov_x;
	static int count_end_mov_y;

	while (1) {

		if (i < NUM_OF_SAMPLES) {
			acc_bias_x += ax;
			acc_bias_y += ay;
			acc_bias_z += az;

			gyro_bias_x += gx;
			gyro_bias_y += gy;
			gyro_bias_z += gz;

			i++;

		} else if (i == NUM_OF_SAMPLES && !offsetCaluclatedAcc) {
			i++;

			acc_bias_x = acc_bias_x / NUM_OF_SAMPLES;   //zbieranie offsetow acc
			acc_bias_y = acc_bias_y / NUM_OF_SAMPLES;
			acc_bias_z = acc_bias_z / NUM_OF_SAMPLES;

			gyro_bias_x = gyro_bias_x / NUM_OF_SAMPLES; //zbieranie offsetow acc
			gyro_bias_y = gyro_bias_y / NUM_OF_SAMPLES;
			gyro_bias_z = gyro_bias_z / NUM_OF_SAMPLES;

			printf("ACC offset calculated :X:%5d  Y:%5d  Z:%5d \n",
					(int) (acc_bias_x * 1000), (int) (acc_bias_y * 1000),
					(int) (acc_bias_z * 1000));
			printf("GYRO offset calculated :X:%5d Y:%5d Z:%5d \n",
					(int) (gyro_bias_x * 1000), (int) (gyro_bias_y * 1000),
					(int) (gyro_bias_z * 1000));
			offsetCaluclatedAcc = true;

		}

		else {
			j++;
			ax_avg += ax;
			ay_avg += ay;
			az_avg += az;  //usrednianie probek

			if (j > 16) {
				j = 0;
				ax_avg = ax_avg / 16;
				ay_avg = ay_avg / 16;
				az_avg = az_avg / 16;

				//KOMPENSACJA WYCHYLENIA, ZAKĹ�UCAJÄ„CEGO ODCZYC Z ACC
				float ax_avg_2 = ax_avg;

				if ((1000 * ax_avg) > 30 && (1000 * ax_avg) < -30) //okno nieczulosci na szumy mechaniczne
						{
					ax_avg = 0;

				}

				if ((1000 * ay_avg) > 30 && (1000 * ay_avg) < -30) {
					ay_avg = 0;

				}

				//calkowanie metoda trapezow

				velocity[0] += (ax_avg + ax_1) * dt * 8 * GRAVITY;  // 0.5 *16
				velocity[1] += (ay_avg + ay_1) * dt * 8 * GRAVITY;
				ax_1 = ax_avg;
				ay_1 = ay_avg;       //  zapamietanie poprzednich wartosci

				//wykyrwanie tak, zawanego konca ruchu

				//os_X
				if (abs(1000 * ax_avg) < 20) {
					count_end_mov_x++;

				} else {
					count_end_mov_x = 0;
				}

				if (count_end_mov_x > 3) {
					velocity[0] = 0;
					velocity_1[0] = 0;

				}

				// os_Y

				if (abs(1000 * ay_avg) < 30) {
					count_end_mov_y++;

				} else {
					count_end_mov_y = 0;
				}

				if (count_end_mov_y > 3) {
					velocity[1] = 0;
					velocity_1[1] = 0;
				}

				position[0] += (velocity[0] + velocity_1[0]) * dt * 8;

				position[1] += (velocity[1] + velocity_1[1]) * dt * 8;

				velocity_1[1] = velocity[1];
				velocity_1[0] = velocity[0];
				ax_avg = 0;
				ay_avg = 0;

				if (velAndPosPrint) {

					printf("VEL X: %3.5f Y: %3.5f Z: %3.3f  %4d %4d m/s \n",
							velocity[0], velocity[1], velocity_1[0],
							(int) (1000 * ax_avg_2),
							(int) (1000 * sin(suv.roll) * az));
					printf("PX:%3.2f m  PY:%3.2f m \n", position[0],
							position[1]);
				}

			}
			velocity[2] += az * dt;

		}

		delay(delay_time);

	}
}

//Watek do wysylania danych z IMU

void *imuDataPrint(void *arg) {
	int delay_time;
	delay_time = (int) arg;
	float dt = delay_time / 1000.0;

	static float pitch = 0;
	static float roll = 0;
	static float yaw = 0;

	bool yawGyroGetFirst = 1;
	static uint32_t time_start = millis();
	while (1) {

		if ((time_start + 15e3 < millis()) && yawGyroGetFirst) { //10s
			yawGyro = yawMag;
			yawGyroGetFirst = 0;
			printf("Gyro yaw settled\n\n");
		}

		///////////// z akcelerometru

		pitchAcc = atan2(ay, (sqrt((ax * ax) + (az * az))));
		rollAcc = atan2(-ax, (sqrt((ay * ay) + (az * az))));
		// z magnetometru

		float Yh = (my * cos(roll)) - (mz * sin(roll));
		float Xh = (mx * cos(pitch)) + (my * sin(roll) * sin(pitch))
				+ (mz * cos(roll) * sin(pitch));

		yawMag = atan2(Yh, Xh);

		//calkowanie

		if (offsetCaluclatedAcc) {

			pitchGyro += (gx * M_PI / 180.0) * dt;
			rollGyro += (gy * M_PI / 180.0) * dt;

			//mechanical filtering
			if (abs(gz) < 2)
				gz = 0;

			yawGyro += (gz * gyroYawFactor * M_PI / 180.0) * dt;

			//if(yawGyro >2*M_PI){ yawGyro= 2*0*M_PI;}
			//if (yawGyro < 0.0) yawGyro += 2.0*M_PI;

			//yawGyro=angle360(yawGyro);

		}

		pitch += (gx * M_PI / 180) * dt;
		roll += (gy * M_PI / 180) * dt;

		yaw += (gz * gyroYawFactor * M_PI / 180.0) * dt;

		float a = 0.04;

		//complementary filter
		pitch = pitchAcc * (a) + pitch * (1 - a);
		roll = rollAcc * (a) + roll * (1 - a);
		yaw = yawMag * (a) + yaw * (1 - a);

		//do obiektu przypsiujemy by miec globalny zasieg
		suv.yaw = yaw;
		suv.roll = roll;
		suv.pitch = pitch;

		eul(0, 0) = pitch;
		eul(1, 0) = roll;
		eul(2, 0) = yaw;

		q2 = Eul2Quat(eul);

		// printf( "P: %4.2f R: %4.2f Y: %4.2f \n",pitchGyro*57.3, rollGyro*57.3,yawGyro*57.3);

		if (eulerPrint) {
			//printf( "P: %4.2f R: %4.2f Y: %4.2f \n",pitch*57.3, roll*57.3,yaw*57.3);
			printf("P: %5.2f R: %5.2f Y: %4.2f GY:%4.2f MY: %4.2f \n",
					pitch * 57.3, roll * 57.3, angle360(yaw) * 57.3,
					angle360(yawGyro) * 57.3, angle360(yawMag) * 57.3);
		}
		//printf(" YG: %4.2f  YC: %4.2f   YM:%4.2f \n",yawGyro*57.3, yaw*57.3, yawMag*57.3);
		//printf( "P: %4d R: %4d Y: %4d \n",(int)(pitchAcc*57.3),(int) (rollAcc*57.3),(int)(yawMag*57.3));
		// Print temperature in degrees Centigrade

		//printf("Gyro temperature is %+1.1f degrees C\n", temperature);
		delay(delay_time);
	}
	return NULL;

}

void *imuGPSPrint(void *arg) {  //GPS 
	int delay_time = (int) arg;
	while (1) {
		while (serialDataAvail(fd)) {
			printf("%c", serialGetchar(fd));
			fflush(stdout);
		}
		delay(delay_time);
	}
	return NULL;
}
void *RadioDataProces(void *arg) {
	int delay_time = (int) arg;
	char data[50];
	int daneRx[8];

	if ((fd_Radio = serialOpen("/dev/ttyUSB0", 115200)) < 0) {
		fprintf(stderr, "Unable to open serial device: %s\n", strerror(errno));

	}

	uint16_t s;
	while (1) {

		while (serialDataAvail(fd_Radio)) {
			read(fd_Radio, &s, 1);
			if (s == 10) {  //nowa linijka
				read(fd_Radio, data, 50);
				sscanf(data, "%d %d %d %d %d %d %d", &daneRx[0], &daneRx[1],
						&daneRx[2], &daneRx[3], &daneRx[4], &daneRx[5],
						&daneRx[6]);
				if (daneRx[0] == 1155) { //ramka startu

					RadioChannel.armed = daneRx[1];
					RadioChannel.forward_backward = daneRx[2];
					RadioChannel.left_right = daneRx[3];
					RadioChannel.auto_mode = daneRx[4];
					RadioChannel.send_gps_pos = daneRx[5];
					RadioChannel.set_point = daneRx[6];

				}

			}

			// printf("%5d %5d %5d %5d %5d %5d %5d  \n",daneRx[0],daneRx[1],daneRx[2],daneRx[3],daneRx[4],daneRx[5],daneRx[6]);
		}
		delay(delay_time);
	}
	return NULL;

}

void *adcRead(void *arg) {
	int delay_time = (int) arg;
	while (1) {
		for (int chan = 0; chan < 8; ++chan) {
			aRead[chan] = analogRead(BASE + chan);
			//printf("%5d",aRead[chan]);
		}
		//printf("\n");
		delay(delay_time);
	}
	return NULL;
}

void *motorsWrite(void *arg) {
	int delay_time = (int) arg;
	uint8_t GYRO = 1;
	float dT = delay_time / 1000.0;
	while (1) {
		if (RadioChannel.armed) {

			suv.speed_motor_A = (RadioChannel.forward_backward - 512)
					+ (RadioChannel.left_right - 512);
			//ogranicznie pred
			if (suv.speed_motor_A > 510)
				suv.speed_motor_A = 510;
			if (suv.speed_motor_A < -510)
				suv.speed_motor_A = -510;

			suv.speed_motor_B = (RadioChannel.forward_backward - 512)
					- (RadioChannel.left_right - 512);
			if (suv.speed_motor_B > 510)
				suv.speed_motor_B = 510;
			if (suv.speed_motor_B < -510)
				suv.speed_motor_B = -510;

			//wystawienie wyjsc pwm

			//silinki 1
			if (suv.speed_motor_A > 30) { //prog zalaczenia

				gpioPWM(MOTOR_A_FORWARD, suv.speed_motor_A * 0.5);
				gpioPWM(MOTOR_A_BACKWARD, 0);
			} else if (suv.speed_motor_A < -30) {

				gpioPWM(MOTOR_A_FORWARD, 0);
				gpioPWM(MOTOR_A_BACKWARD, -(suv.speed_motor_A * 0.5));
			} else {
				gpioPWM(MOTOR_A_FORWARD, 0);
				gpioPWM(MOTOR_A_BACKWARD, 0);
			}

			//silnik 2
			if (suv.speed_motor_B > 30) { //prog zalaczenia

				gpioPWM(MOTOR_B_FORWARD, suv.speed_motor_B * 0.5);
				gpioPWM(MOTOR_B_BACKWARD, 0);
			} else if (suv.speed_motor_B < -30) {

				gpioPWM(MOTOR_B_FORWARD, 0);
				gpioPWM(MOTOR_B_BACKWARD, -(suv.speed_motor_B * 0.5));
			} else {
				gpioPWM(MOTOR_B_FORWARD, 0);
				gpioPWM(MOTOR_B_BACKWARD, 0);
			}

			//printf(" %4d  %4d  \n", suv.speed_motor_A,suv.speed_motor_B);

		} else if (GYRO) {
			yaw_dest = 90;
			float yaw_e = yaw_dest - (int) (suv.yaw * 57.3);
			float Kp = 0.7;
			float Ki = 0.5;

			float correction = reg_PI(yaw_e, Kp, Ki, 500, -500, dT);

			suv.speed_motor_A = correction;
			suv.speed_motor_B = -correction;

		}

		else {

			gpioPWM(MOTOR_B_FORWARD, 0);
			gpioPWM(MOTOR_B_BACKWARD, 0);
			gpioPWM(MOTOR_A_FORWARD, 0);
			gpioPWM(MOTOR_A_BACKWARD, 0);

		}
		if (motorsPrint)
			printf(" %4d  %4d  %4f  \n", suv.speed_motor_A, suv.speed_motor_B,
					suv.yaw * 57.3);

		delay(delay_time);
	}
	return NULL;
}

#define RAD_TO_DEG  180.0/M_PI
#define DEG_TO_RAD  M_PI/180.0

void *logDataToFile(void *arg) {

	int delay_time = (int) arg;
	FILE * logFile = fopen("dataLog.txt", "w");
	FILE * logFileRaw = fopen("dataRawLog.txt", "w");

	while (1) {

		fprintf(logFile, "%8d %4.3f %4.3f %4.3f ", millis(),
				suv.roll * RAD_TO_DEG, suv.pitch * RAD_TO_DEG,
				suv.yaw * RAD_TO_DEG); //roll. pitch ,yaw, degrees from complemetnary
		//fprintf(logFile, "%4.3f %4.3f %4.3f", suv.roll*RAD_TO_DEG, suv.pitch*RAD_TO_DEG, suv.yaw*RAD_TO_DEG);

		fprintf(logFile, " %4.3f %4.3f  ", rollGyro * RAD_TO_DEG,
				rollAcc * RAD_TO_DEG);
		fprintf(logFile, " %4.3f %4.3f  ", pitchGyro * RAD_TO_DEG,
				pitchAcc * RAD_TO_DEG);
		fprintf(logFile, " %4.3f %4.3f  ", yawGyro * RAD_TO_DEG,
				yawMag * RAD_TO_DEG);
		fprintf(logFile, " \n\r");

		// "%8d %4 %4.3f %4.3f",millis(), suv.roll*RAD_TO_DEG, suv.pitch*RAD_TO_DEG, suv.yaw*RAD_TO_DEG);
		fprintf(logFileRaw, "%8d %5d  %5d  %5d   ", millis(), (int) (1000 * ax),
				(int) (1000 * ay), (int) (1000 * az));   //mg
		fprintf(logFileRaw, "%+2.2f  %+2.2f   %+2.2f ", gx, gy, gz);  //deg/s
		fprintf(logFileRaw, "%4d  %4d  %4d \r\n", (int) mx, (int) my, (int) mz); //mG

		//printf("\n");
		delay(delay_time);
	}
	return NULL;

}

float constrainAngle180(float dta) {
	if (dta > M_PI)
		dta -= (M_PI * 2.0);
	if (dta < -M_PI)
		dta += (M_PI * 2.0);
	return dta;
}

void *imuAHRS_EKF(void *arg) {
	bool firstTime = true;

	int delay_time = (int) arg;
	filtrEKF.setInitializationDuration(6e6); //60s

	while (1) {

		if (filtrEKF.update(gx * DEG_TO_RAD, gy * DEG_TO_RAD, gz * DEG_TO_RAD,
				ax, ay, az, mx / 1000.0, my / 1000.0, mz / 1000.0)) {
			if (firstTime) {

				printf("EKF init finished! \n");
				firstTime = false;
			}

			if (printEKF) {
				ekf.pitch = filtrEKF.getPitch_rad() * RAD_TO_DEG;
				ekf.roll = constrainAngle180(
						-M_PI + filtrEKF.getRoll_rad())*RAD_TO_DEG;
				ekf.yaw = 360 - filtrEKF.getHeading_rad() * RAD_TO_DEG;

				printf("%5.3f %5.3f %5.3f  %5.3f \n", ekf.pitch, ekf.roll,
						ekf.yaw, angle360(suv.yaw) * RAD_TO_DEG); //chyba

			}
		}
		delay(delay_time);
	}
	return NULL;

}

void setup() {
	////////////////////////////SETUP///////////////////////////////////

	// Setup WirinPi
	wiringPiSetup();
	pigpioInit();
	//inicjalizacja GPS
	gpsInit();
	// inicjalizacja przetwornika ADC
	mcp3004Setup(BASE, SPI_CHAN);
	// ustawienie pinow PWM do sterowania silnikami
	gpioSetMode(MOTOR_A_BACKWARD, PI_OUTPUT);
	gpioSetMode(MOTOR_A_FORWARD, PI_OUTPUT);

	gpioSetMode(MOTOR_B_BACKWARD, PI_OUTPUT);
	gpioSetMode(MOTOR_B_FORWARD, PI_OUTPUT);
	// Start the MPU9250
	switch (imu.begin()) {

	case MPUIMU::ERROR_IMU_ID:
		// errmsg("Bad IMU device ID");
		break;
	case MPUIMU::ERROR_MAG_ID:
		//errmsg("Bad magnetometer device ID");
		break;
	case MPUIMU::ERROR_SELFTEST:
		//errmsg("Failed self-test");
		break;
	default:
		printf("MPU9250 online!\n");
		break;
	}


	//Tworzenie wątków do wykonwyania poszczególnych zadań. Funkcja pthread_create przyjmuje jako argument wsakźnik
	//do funkcji oraz ilosc ms co jaką wykonywane jest dane zadanie

	if (pthread_create(&watekRadioData, NULL, RadioDataProces, (void*) 50)) { //20Hz
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekImuEuler, NULL, imuDataPrint, (void*) 20)) { //50Hz
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekImuReadData, NULL, imuRead, (void*) 4)) { //330Hz
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekGpsReadData, NULL, imuGPSPrint, (void*) 10)) { //100Hz
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekMotorsWrite, NULL, motorsWrite, (void*) 20)) { //50Hz
		printf("Błąd przy tworzeniu watku");
		abort();
	}
	if (pthread_create(&watekADCReadData, NULL, adcRead, (void*) 100)) {
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekPostion, NULL, imuPositionCalc, (void*) 5)) {
		printf("Błąd przy tworzeniu watku");
		abort();
	}

	if (pthread_create(&watekLogData, NULL, logDataToFile, (void*) 20)) {
		printf("Błąd przy tworzeniu watku");
		abort();
	}



	wiringPiISR(intPin, INT_EDGE_RISING, &myinthandler); // define interrupt for intPin output of MPU9250
}
//////////////////////////////////////////////PETLA GLOWNA////////////////
void loop()   // menu tekstowe w pętli głównej do obsługi poszczególnych funkcji kontolno kalibracyjnych
{

	printf("Podaj czynnosc :\n");
	printf(" [1] Podawanie kÄ…ta YAW\n");
	printf(" [2] RÄ™czne nastawianie wartoĹ›ci predkosci silnikow\n");
	printf(" [3] Konfiguracja regulatorw PI\n");
	printf(" [4] Kalibracja magnetometru\n");
	printf(" [5] Katy Eulera -filtr komplemetnarny \n");
	printf(" [6] Polozenie X, Y, Z\n");
	printf(" [7] Motors values <0 - 512> \n");
	printf(" [8] Dane RAW\n");
	printf(" [9] Katy Eulera -EKF\n");
	printf(" [10] Wartosci kalibracyjne z pliku\n");

	//oczekiwanie na polecenie
	scanf("%d", &polecenie);

	switch (polecenie) {

	case 1:
		int polecenieYaw;
		printf("Podaj wartosc kata YAW ( -90 do 90)\n");
		scanf("%d", &polecenieYaw);
		printf("Podales kat yaw : %d \n", polecenieYaw);
		yaw_dest = polecenieYaw;

		break;

	case 0:
		motorsPrint = 0;
		eulerPrint = 0;
		velAndPosPrint = 0;
		senosorsRawPrint = 0;
		printEKF = 0;
		break;

	case 4:
		printf("Rozpoczeto kalibracje magnetometru.... - > wykonuj ruchy osemkowe!\n");
		imu.calibrateMagnetometer();

		yawGyro = suv.yaw;
		imu.readMagBiasAndScale(bx, by, bz, sx, sy, sz);
		printf("Magnetometetr skalibrowany.\n");
		printf("Biases X: %f Y:%f  Z: %f \n", bx, by, bz);
		printf("Scales X: %f Y:%f  Z: %f \n", sx, sy, sz);

		break;

	case 10:
		printf("Mag auto -calibreted\n");
		imu.setMagBiasAndScale(bx, by, bz, sx, sy, sz);
		imu.readMagBiasAndScale(bx, by, bz, sx, sy, sz);
		printf("Magnetometetr skalibrowany.\n");
		printf("Biases X: %f Y:%f  Z: %f \n", bx, by, bz);
		printf("Scales X: %f Y:%f  Z: %f \n", sx, sy, sz);
		break;

	case 5:
		eulerPrint = 1;
		break;

	case 6:
		velAndPosPrint = 1;
		break;
	case 7:
		motorsPrint = 1;
		break;
	case 8:
		senosorsRawPrint = 1;
		break;
	case 9:
		printEKF = 1;
		break;
	case 3:

		if (pthread_create(&watekEKF, NULL, imuAHRS_EKF, (void*) 10)) {
			printf("Blad przy tworzeniu watku");
			abort();
		}
		break;
	case 33:
		calibrateGyroyYaw();
		break;

	default: /
		break;

	}

}

