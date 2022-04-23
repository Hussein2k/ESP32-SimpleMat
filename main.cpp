#include <Arduino.h>
#include "templateMatrix.h"
#include <math.h>
#include <esp_adc_cal.h>

#define M_PI 3.141592653
#define ADC_PIN 34



float fast_rand(float a, float b) {
  return a + (rand() / float(RAND_MAX)) * (b - a);
}

void setup() {

#ifdef SERIAL_DBG
  Serial.begin(9600);
#endif
  pinMode(ADC_PIN, INPUT);

}
inline float read_val() {
  return analogRead(ADC_PIN) * 0.05f;
}
void loop() {

  float dt = 0.001;
  float w = 2 * M_PI;
  const int samples = 10000;
  srand(time(0));

  auto A = Matrix<float, 2, 2>({
       {1,dt},
       {-dt * w * w,1}
    });
  auto H = Matrix<float, 2, 2>({
      {1,0},
      {0,0}
    });

  auto Q = Matrix<float, 2, 2>({
      {1e-3,1e-4},
      {1e-4,1e-3}
    });
  auto P = Matrix<float, 2, 2>({
      {1,1},
      {1,1}
    });
  auto R = Matrix<float, 2, 2>({
      {0.8,0.64},
      {0.64,0.8}
    });


  auto I = P;
  Identity(I);

  //Priori
  auto X = Matrix<float, 2, 1>({
      {0},{w}
    });
  auto X_ = A * X;
  auto P_ = A * P * Transpose(A) + Q;

  //posteriori
  auto Ht = Transpose(H);
  auto K = P_ * Ht * Inverse(H * P_ * Ht + R);
  auto M = Matrix<float, 2, 1>({ {read_val()},{0} });//Measurements

  X = X_ + K * (M - H * X_);
  P = (I - K * H) * P_;


  for (;;) {
    //Priori
    X_ = A * X;
    P_ = A * P * Transpose(A) + Q;

    //posteriori
    Ht = Transpose(H);
    K = P_ * Ht * Inverse(H * P_ * Ht + R);
    M = Matrix<float, 2, 1>({ {read_val()},{0} });//Measurements
    X = X_ + K * (M - H * X_);
    P = (I - K * H) * P_;
#ifdef SERIAL_DBG
    Serial.println(*X.get_elem(0, 0));
#endif
  }


}