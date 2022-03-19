#ifndef PID_H
#define PID_H


class PID
{
public:
    float Ki; // 积分
    float Kp; // 比例
    float Kf; // 微分

    float PreviousErr = 0.0;
    float CurrentErr  = 0.0;
    float DiffErr     = 0.0;
    float InteErr     = 0.0;
    float Output      = 0.0;

    PID();
    PID(float Ki, float Kp, float Kf);
    void CalculateError(float Input); // Must Call it at
    void CalculateDiff(); // 计算微分量
    void CalculateInte(); // 计算积分量
    void Update(float Input);
};



#endif // PID_H
