#include "PID.h"

PID::PID()
{

}

PID::PID(float Ki, float Kp, float Kf){
    this->Ki = Ki;
    this->Kp = Kp;
    this->Kf = Kf;
}

void PID::CalculateError(float Input) // Must Call it at 500Hz
{
    this->PreviousErr = this->CurrentErr;
    this->CurrentErr = Input - Output;
}

void PID::CalculateDiff() // 计算微分量
{
    this->DiffErr = this->CurrentErr - this->PreviousErr;
}

void PID::CalculateInte() // 计算积分量
{
    this->InteErr += this->CurrentErr;
}

void PID::Update(float Input){
    CalculateError(Input);
    CalculateDiff();
    CalculateInte();
    Output += (Ki * this->InteErr + Kp * this->CurrentErr + Kf * this->DiffErr);
}
