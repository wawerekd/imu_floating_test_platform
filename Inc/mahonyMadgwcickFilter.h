
// do zadeklarowania w main.c

extern float q[4];
extern  float beta;
extern  float deltat;
extern  float Kp;
extern  float Ki;
extern  float eInt[4];

//funckcje

void MadgwickQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
void MahonyQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);


