## IKFoM 
**IKFoM** （在流形上迭代卡尔曼滤波器）是一个计算高效且方便的工具包，旨在为各种机器人系统部署迭代卡尔曼滤波器，特别是那些在高维流形上运行的系统。它实现了一种流形嵌入的卡尔曼滤波器，可以将流形结构与系统描述分离，并且用户只需通过定义系统的规范形式并相应地调用各个步骤即可使用。    
当前的实现支持在流形 <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{R}^m\times&space;SO(3)\times\cdots\times&space;SO(3)\times\mathbb{S}^2\times\cdots\times\mathbb{S}^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{R}^m\times&space;SO(3)\times\cdots\times&space;SO(3)\times\mathbb{S}^2\times\cdots\times\mathbb{S}^2" title="\mathbb{R}^m\times SO(3)\times\cdots\times SO(3)\times\mathbb{S}^2\times\cdots\times\mathbb{S}^2" /></a> 上进行完整的迭代卡尔曼滤波，及其任意子流形，并能够在必要时扩展到其他类型的流形。


**贡献者**

[Dongjiao He](https://github.com/Joanna-HE)

**演示视频**: https://youtu.be/sz_ZlDkl6fA

## 1. 先决条件

### 1.1. **Eigen && Boost**
Eigen  >= 3.3.4, 参考 [Eigen Installation](http://eigen.tuxfamily.org/index.php?title=Main_Page).

Boost >= 1.65.

## 2. 当测量值为恒定尺寸和类型时使用
克隆储存库:

```
    git clone https://github.com/hku-mars/IKFoM.git
```

1. 包含必要的头文件:
```
#include<esekfom/esekfom.hpp>
```
2. 选择并实例化原始流形:
```
    typedef MTK::SO3<double> SO3; // scalar type of variable: double
    typedef MTK::vect<3, double> vect3; // dimension of the defined Euclidean variable: 3
    typedef MTK::S2<double, 98, 10, 1> S2; // length of the S2 variable: 98/10; choose e1 as the original point of rotation: 1
```
3. 将系统状态、输入和测量构建为由原始流形组成的复合流形:
``` 
MTK_BUILD_MANIFOLD(state, // name of compound manifold: state
((vect3, pos)) // ((primitive manifold type, name of variable))
((vect3, vel))
((SO3, rot))
((vect3, bg))
((vect3, ba))
((S2, grav))
((SO3, offset_R_L_I))
((vect3, offset_T_L_I)) 
);
```
4. 实现向量场 <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" title="\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)" /></a> 的定义 <a href="https://latex.codecogs.com/svg.image?\mathbf{x}_{k&plus;1}&space;=&space;\mathbf{x}_k\oplus\Delta&space;t\mathbf{f}(\mathbf{x}_k,&space;\mathbf{u}_k,&space;\mathbf{w}_k);\hat{\mathbf{x}}_{k&plus;1}&space;=&space;\hat{\mathbf{x}}_k\oplus\Delta&space;t\mathbf{f}(\hat{\mathbf{x}}_k,&space;\mathbf{u}_k,&space;\mathbf{0})"><see here>, 及其微分 <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{u}, \mathbf{0}\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)}{\partial\mathbf{w}}" /></a>, 其中 w=0 可以省略:
```
Eigen::Matrix<double, state_length, 1> f(state &s, const input &i) {
	Eigen::Matrix<double, state_length, 1> res = Eigen::Matrix<double, state_length, 1>::Zero();
	res(0) = s.vel[0];
	res(1) = s.vel[1];
	res(2) = s.vel[2];
	return res;
}
Eigen::Matrix<double, state_length, state_dof> df_dx(state &s, const input &i) //notice S2 has length of 3 and dimension of 2 {
	Eigen::Matrix<double, state_length, state_dof> cov = Eigen::Matrix<double, state_length, state_dof>::Zero();
	cov.template block<3, 3>(0, 12) = Eigen::Matrix3d::Identity();
	return cov;
}
Eigen::Matrix<double, state_length, process_noise_dof> df_dw(state &s, const input &i) {
	Eigen::Matrix<double, state_length, process_noise_dof> cov = Eigen::Matrix<double, state_length, process_noise_dof>::Zero();
	cov.template block<3, 3>(12, 3) = -s.rot.toRotationMatrix();
	return cov;
}
```
这些函数将在 ekf state predict 期间被调用

5. 实现输出方程 <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> 及其微分 <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
measurement h(state &s, bool &valid) // the iteration stops before convergence whenever the user set valid as false
{
	if (condition){ valid = false; 
	} // other conditions could be used to stop the ekf update iteration before convergence, otherwise the iteration will not stop until the condition of convergence is satisfied.
	measurement h_;
	h_.position = s.pos;
	return h_;
}
Eigen::Matrix<double, measurement_dof, state_dof> dh_dx(state &s, bool &valid) {} 
Eigen::Matrix<double, measurement_dof, measurement_noise_dof> dh_dv(state &s, bool &valid) {}
```
这些函数将在 ekf 状态更新期间被调用

6. 实例化 **esekf** 对象 **kf** 并使用初始或默认状态和协变对其进行初始化。

(1) 初始化状态和协方差:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof>::cov init_P;
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof> kf(init_state,init_P);
```
(2) 默认状态和协方差
```
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof> kf;
```
其中 **process_noise_dof** 是进程噪声的维度, 类型为 std int, 因此 **measurement_noise_dof**.

7. 将定义的模型、std int 最大迭代次数 **Maximum_iter**, 和用于测试收敛 **epsi** 的 std 数组交付到 **esekf** 对象中:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init(f, df_dx, df_dw, h, dh_dx, dh_dv, Maximum_iter, epsi);
```
8. 在运行时中，一旦在最后一次传播或更新后收到 **in** 或测量 **z** 后 **dt** 就会执行传播:
```
kf.predict(dt, Q, in); // process noise covariance: Q, an Eigen matrix
```
9. 收到测量 **z** 后, 将执行迭代更新:
```
kf.update_iterated(z, R); // measurement noise covariance: R, an Eigen matrix
```
*备注(1):*
- 我们还将输出方程及其微分组合成一个联合函数，其用法与上述步骤 1-4 相同，步骤 5-9 如下所示.

5. 实现输出方程 <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> 及其微分 <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
measurement h_share(state &s, esekfom::share_datastruct<state, measurement, measurement_noise_dof> &share_data) 
{
	if(share_data.converge) {} // this value is true means iteration is converged 
	if(condition) share_data.valid = false; // the iteration stops before convergence when this value is false if other conditions are satified
	share_data.h_x = H_x; // H_x is the result matrix of the first differentiation 
	share_data.h_v = H_v; // H_v is the result matrix of the second differentiation
	share_data.R = R; // R is the measurement noise covariance
	share_data.z = z; // z is the obtained measurement 

	measurement h_;
	h_.position = s.pos;
	return h_;
}
```
这个函数将在 ekf state 更新期间调用，输出函数及其导数、测量和测量噪声将从这个联合函数中获得。

6. 实例化 **esekf** 对象 **kf** 并使用初始或默认状态和协方差对其进行初始化.

(1) 初始状态和协方差:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof>::cov init_P;
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof> kf(init_state,init_P);
```
(2) 默认状态和协方差:
```
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof> kf;
```
1. Deliver the defined models, std int maximum iteration numbers **Maximum_iter**, and the std array for testing convergence **epsi** into the **esekf** object:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init_share(f, df_dx, df_dw, h_share, Maximum_iter, epsi);
```
1. In the running time, once an input **in** or a measurement **z** is received **dt** after the last propagation or update, a propagation is executed:
```
kf.predict(dt, Q, in); // process noise covariance: Q
```
1. Once a measurement **z** is received, an iterated update is executed:
```
kf.update_iterated_share();
```

*Remarks(2):*
- The value of the state **x** and the covariance **P** are able to be changed by functions **change_x()** and **change_P()**:
```
state set_x;
kf.change_x(set_x);
esekfom::esekf<state, process_noise_dof, input, measurement, measurement_noise_dof>::cov set_P;
kf.change_P(set_P);
```

## 3. Usage when the measurement is an Eigen vector of changing dimension.

Clone the repository:

```
    git clone https://github.com/hku-mars/IKFoM.git
```

1. include the necessary head file:
```
#include<esekfom/esekfom.hpp>
```
2. Select and instantiate the primitive manifolds:
```
    typedef MTK::SO3<double> SO3; // scalar type of variable: double
    typedef MTK::vect<3, double> vect3; // dimension of the defined Euclidean variable: 3
    typedef MTK::S2<double, 98, 10, 1> S2; // length of the S2 variable: 98/10; choose e1 as the original point of rotation: 1
```
3. Build system state and input as compound manifolds which are composed of the primitive manifolds:
``` 
MTK_BUILD_MANIFOLD(state, // name of compound manifold: state
((vect3, pos)) // ((primitive manifold type, name of variable))
((vect3, vel))
((SO3, rot))
((vect3, bg))
((vect3, ba))
((S2, grav))
((SO3, offset_R_L_I))
((vect3, offset_T_L_I)) 
);
```
4. Implement the vector field <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" title="\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)" /></a> that is defined as <a href="https://latex.codecogs.com/svg.image?\mathbf{x}_{k&plus;1}&space;=&space;\mathbf{x}_k\oplus\Delta&space;t\mathbf{f}(\mathbf{x}_k,&space;\mathbf{u}_k,&space;\mathbf{w}_k);\hat{\mathbf{x}}_{k&plus;1}&space;=&space;\hat{\mathbf{x}}_k\oplus\Delta&space;t\mathbf{f}(\hat{\mathbf{x}}_k,&space;\mathbf{u}_k,&space;\mathbf{0})"> <see here>, and its differentiation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{u}, \mathbf{0}\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)}{\partial\mathbf{w}}" /></a>, where w=0 could be left out:
```
Eigen::Matrix<double, state_length, 1> f(state &s, const input &i) {
	Eigen::Matrix<double, state_length, 1> res = Eigen::Matrix<double, state_length, 1>::Zero();
	res(0) = s.vel[0];
	res(1) = s.vel[1];
	res(2) = s.vel[2];
	return res;
}
Eigen::Matrix<double, state_length, state_dof> df_dx(state &s, const input &i) //notice S2 has length of 3 and dimension of 2 {
	Eigen::Matrix<double, state_length, state_dof> cov = Eigen::Matrix<double, state_length, state_dof>::Zero();
	cov.template block<3, 3>(0, 12) = Eigen::Matrix3d::Identity();
	return cov;
}
Eigen::Matrix<double, state_length, process_noise_dof> df_dw(state &s, const input &i) {
	Eigen::Matrix<double, state_length, process_noise_dof> cov = Eigen::Matrix<double, state_length, process_noise_dof>::Zero();
	cov.template block<3, 3>(12, 3) = -s.rot.toRotationMatrix();
	return cov;
}
```
Those functions would be called during ekf state predict

5. Implement the output equation <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> and its differentiation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
Eigen::Matrix<double, Eigen::Dynamic, 1> h(state &s, bool &valid) //the iteration stops before convergence when valid is false {
	if (condition){ valid = false; 
	} // other conditions could be used to stop the ekf update iteration before convergence, otherwise the iteration will not stop until the condition of convergence is satisfied.
	Eigen::Matrix<double, Eigen::Dynamic, 1> h_;
	h_(0) = s.pos[0];
	return h_;
}
Eigen::Matrix<double, Eigen::Dynamic, state_dof> dh_dx(state &s, bool &valid) {} 
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dh_dv(state &s, bool &valid) {}
```
Those functions would be called during ekf state update

6. Instantiate an **esekf** object **kf** and initialize it with initial or default state and covariance.

(1) initial state and covariance:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input>::cov init_P;
esekfom::esekf<state, process_noise_dof, input> kf(init_state,init_P);
```
(2) default state and covariance:
```
esekfom::esekf<state, process_noise_dof, input> kf;
```
where **process_noise_dof** is the dimension of process noise, with the type of std int, and so for **measurement_noise_dof**

7. Deliver the defined models, std int maximum iteration numbers **Maximum_iter**, and the std array for testing convergence **epsi** into the **esekf** object:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init_dyn(f, df_dx, df_dw, h, dh_dx, dh_dv, Maximum_iter, epsi);
```
8. In the running time, once an input **in** or an measurement **z** is received **dt** after the last propagation or update, a propagation is executed:
```
kf.predict(dt, Q, in); // process noise covariance: Q, an Eigen matrix
```
9. Once a measurement **z** is received, an iterated update is executed:
```
kf.update_iterated_dyn(z, R); // measurement noise covariance: R, an Eigen matrix
```
*Remarks(1):*
- We also combine the output equation and its differentiation into an union function, whose usage is the same as the above steps 1-4, and steps 5-9 are shown as follows.
5. Implement the output equation <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> and its differentiation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
Eigen::Matrix<double, Eigen::Dynamic, 1> h_dyn_share(state &s, esekfom::dyn_share_datastruct<double> &dyn_share_data) 
{
	if(dyn_share_data.converge) {} // this value is true means iteration is converged 
	if(condition) share_data.valid = false; // the iteration stops before convergence when this value is false if other conditions are satified
	dyn_share_data.h_x = H_x; // H_x is the result matrix of the first differentiation
	dyn_share_data.h_v = H_v; // H_v is the result matrix of the second differentiation 
	dyn_share_data.R = R; // R is the measurement noise covariance
	dyn_share_data.z = z; // z is the obtained measurement 

	Eigen::Matrix<double, Eigen::Dynamic, 1> h_;
	h_(0) = s.pos[0];
	return h_;
}
This function would be called during ekf state update, and the output function and its derivatives, the measurement and the measurement noise would be obtained from this one union function
```
6. Instantiate an **esekf** object **kf** and initialize it with initial or default state and covariance.
(1) initial state and covariance:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input>::cov init_P;
esekfom::esekf<state, process_noise_dof, input> kf(init_state,init_P);
```
(2) default state and covariance:
```
esekfom::esekf<state, process_noise_dof, input> kf;
```
7. Deliver the defined models, std int maximum iteration numbers **Maximum_iter**, and the std array for testing convergence **epsi** into the **esekf** object:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init_dyn_share(f, df_dx, df_dw, h_dyn_share, Maximum_iter, epsi);
```
8. In the running time, once an input **in** or a measurement **z** is received **dt** after the last propagation or update, a propagation is executed:
```
kf.predict(dt, Q, in); // process noise covariance: Q, an Eigen matrix
```
9. Once a measurement **z** is received, an iterated update is executed:
```
kf.update_iterated_dyn_share();
```

*Remarks(2):*
- The value of the state **x** and the covariance **P** are able to be changed by functions **change_x()** and **change_P()**:
```
state set_x;
kf.change_x(set_x);
esekfom::esekf<state, process_noise_dof, input>::cov set_P;
kf.change_P(set_P);
```

## 4. Usage when the measurement is a changing manifold during the run time.

Clone the repository:

```
    git clone https://github.com/hku-mars/IKFoM.git
```

1. include the necessary head file:
```
#include<esekfom/esekfom.hpp>
```
2. Select and instantiate the primitive manifolds:
```
    typedef MTK::SO3<double> SO3; // scalar type of variable: double
    typedef MTK::vect<3, double> vect3; // dimension of the defined Euclidean variable: 3
    typedef MTK::S2<double, 98, 10, 1> S2; // length of the S2 variable: 98/10; choose e1 as the original point of rotation: 1
```
3. Build system state and input as compound manifolds which are composed of the primitive manifolds:
``` 
MTK_BUILD_MANIFOLD(state, // name of compound manifold: state
((vect3, pos)) // ((primitive manifold type, name of variable))
((vect3, vel))
((SO3, rot))
((vect3, bg))
((vect3, ba))
((S2, grav))
((SO3, offset_R_L_I))
((vect3, offset_T_L_I)) 
);
```
4. Implement the vector field <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)" title="\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)" /></a> that is defined as <a href="https://latex.codecogs.com/svg.image?\mathbf{x}_{k&plus;1}&space;=&space;\mathbf{x}_k\oplus\Delta&space;t\mathbf{f}(\mathbf{x}_k,&space;\mathbf{u}_k,&space;\mathbf{w}_k);\hat{\mathbf{x}}_{k&plus;1}&space;=&space;\hat{\mathbf{x}}_k\oplus\Delta&space;t\mathbf{f}(\hat{\mathbf{x}}_k,&space;\mathbf{u}_k,&space;\mathbf{0})"> <see here>, and its differentiation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{u},&space;\mathbf{0}\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{u}, \mathbf{0}\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\mathbf{f}\left(\mathbf{x},&space;\mathbf{u},&space;\mathbf{w}\right)}{\partial\mathbf{w}}" title="\frac{\partial\mathbf{f}\left(\mathbf{x}, \mathbf{u}, \mathbf{w}\right)}{\partial\mathbf{w}}" /></a>, where w=0 could be left out:
```
Eigen::Matrix<double, state_length, 1> f(state &s, const input &i) {
	Eigen::Matrix<double, state_length, 1> res = Eigen::Matrix<double, state_length, 1>::Zero();
	res(0) = s.vel[0];
	res(1) = s.vel[1];
	res(2) = s.vel[2];
	return res;
}
Eigen::Matrix<double, state_length, state_dof> df_dx(state &s, const input &i) //notice S2 has length of 3 and dimension of 2 {
	Eigen::Matrix<double, state_length, state_dof> cov = Eigen::Matrix<double, state_length, state_dof>::Zero();
	cov.template block<3, 3>(0, 12) = Eigen::Matrix3d::Identity();
	return cov;
}
Eigen::Matrix<double, state_length, process_noise_dof> df_dw(state &s, const input &i) {
	Eigen::Matrix<double, state_length, process_noise_dof> cov = Eigen::Matrix<double, state_length, process_noise_dof>::Zero();
	cov.template block<3, 3>(12, 3) = -s.rot.toRotationMatrix();
	return cov;
}
```
Those functions would be called during ekf state predict

5. Implement the differentiation of the output equation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
Eigen::Matrix<double, Eigen::Dynamic, state_dof> dh_dx(state &s, bool &valid) {} //the iteration stops before convergence when valid is false
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dh_dv(state &s, bool &valid) {}
```
Those functions would be called during ekf state update

6. Instantiate an **esekf** object **kf** and initialize it with initial or default state and covariance.

(1) initial state and covariance:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input>::cov init_P;
esekfom::esekf<state, process_noise_dof, input> kf(init_state,init_P);
```
(2)
```
esekfom::esekf<state, process_noise_dof, input> kf;
```
Where **process_noise_dof** is the dimension of process noise, of type of std int

7. Deliver the defined models, std int maximum iteration numbers **Maximum_iter**, and the std array for testing convergence **epsi** into the **esekf** object:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init_dyn_runtime(f, df_dx, df_dw, dh_dx, dh_dv, Maximum_iter, epsi);
```
8. In the running time, once an input **in** or a measurement **z** is received **dt** after the last propagation or update, a propagation is executed:
```
kf.predict(dt, Q, in); // process noise covariance: Q
```
9. Once a measurement **z** is received, build system measurement as compound manifolds following step 3 and implement the output equation <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> :
```
measurement h(state &s, bool &valid) //the iteration stops before convergence when valid is false
{
	if (condition) valid = false; // the update iteration could be stopped when the condition other than convergence is satisfied
	measurement h_;
	h_.pos = s.pos;
	return h_;
}
``` 
then an iterated update is executed:
```
kf.update_iterated_dyn_runtime(z, R, h); // measurement noise covariance: R, an Eigen matrix
```
*Remarks(1):*
- We also combine the output equation and its differentiation into an union function, whose usage is the same as the above steps 1-4, and steps 5-9 are shown as follows.
5. Instantiate an **esekf** object **kf** and initialize it with initial or default state and covariance.

(1) initial state and covariance:
```
state init_state;
esekfom::esekf<state, process_noise_dof, input>::cov init_P;
esekfom::esekf<state, process_noise_dof, input> kf(init_state,init_P);
```
(2) default state and covariance:
```
esekfom::esekf<state, process_noise_dof, input> kf;
```
6. Deliver the defined models, std int maximum iteration numbers **Maximum_iter**, and the std array for testing convergence **epsi** into the **esekf** object:
```
double epsi[state_dof] = {0.001};
fill(epsi, epsi+state_dof, 0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is converged
kf.init_dyn_runtime_share(f, df_dx, df_dw, Maximum_iter, epsi);
```
7. In the running time, once an input **in** or a measurement **z** is received **dt** after the last propagation or update, a propagation is executed:
```
kf.predict(dt, Q, in); // process noise covariance: Q. an Eigen matrix
```
8. Once a measurement **z** is received, build system measurement as compound manifolds following step 3 and implement the output equation <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)" title="\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)" /></a> and its differentiation <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x},&space;\mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}\boxplus\delta\mathbf{x}, \mathbf{0}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\delta\mathbf{x}}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial\left(\mathbf{h}\left(\mathbf{x},&space;\mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" title="\frac{\partial\left(\mathbf{h}\left(\mathbf{x}, \mathbf{v}\right)\boxminus\mathbf{h}\left(\mathbf{x},\mathbf{0}\right)\right)}{\partial\mathbf{v}}" /></a>:
```
measurement h_dyn_runtime_share(state &s, esekfom::dyn_runtime_share_datastruct<double> &dyn_runtime_share_data) 
{
	if(dyn_runtime_share_data.converge) {} // this value is true means iteration is converged 
	if(condition) dyn_runtime_share_data.valid = false; // the iteration stops before convergence when this value is false, if conditions other than convergence is satisfied
	dyn_runtime_share_data.h_x = H_x; // H_x is the result matrix of the first differentiation 
	dyn_runtime_share_data.h_v = H_v; // H_v is the result matrix of the second differentiation
	dyn_runtime_share_data.R = R; // R is the measurement noise covariance
	
	measurement h_;
	h_.pos = s.pos;
	return h_;
}
```
This function would be called during ekf state update, and the output function and its derivatives, the measurement and the measurement noise would be obtained from this one union function

then an iterated update is executed:
```
kf.update_iterated_dyn_runtime_share(z, h_dyn_runtime_share);
```

*Remarks(2):*
- The value of the state **x** and the covariance **P** are able to be changed by functions **change_x()** and **change_P()**:
```
state set_x;
kf.change_x(set_x);
esekfom::esekf<state, process_noise_dof, input>::cov set_P;
kf.change_P(set_P);
```

## 5. Run samples
Clone the repository:

```
    git clone https://github.com/hku-mars/IKFoM.git
```
Sample 1. In the **Sample_1** file folder, there is the scource code that applys the **IKFoM** on the original source code from [FAST LIO](https://github.com/hku-mars/FAST_LIO). Please follow the README.md shown in that repository excepting the step **2. Build**, which is modified as:
```
cd ~/catkin_ws/src
cp -r ~/IKFoM/Samples/FAST_LIO FAST_LIO
cd ..
catkin_make
source devel/setup.bash
```

Livox Avia rosbag: Can be downloaded from [here](https://connecthkuhk-my.sharepoint.com/:f:/g/personal/hdj65822_connect_hku_hk/Eu1wAc34q1ZDvHCwe2ajAK4BsH82-c3yTldI9Blp6_rivg?e=41sEEw)

For the indoor bag, run:
```
roslaunch fast_lio mapping_avia.launch
rosbag play YOUR_DOWNLOADED.bag
```

For the outdoor bag, run:
```
roslaunch fast_lio mapping_avia_outdoor.launch
rosbag play YOUR_DOWNLOADED.bag
```

Sample 2. In the **Sample_2** file folder, there is the scource code that applys the **IKFoM** on the original source code from [LINS](https://github.com/ChaoqinRobotics/LINS---LiDAR-inertial-SLAM). Please follow the README.md shown in that repository excepting the step **Compile**, which is modified as:
```
cd ~/catkin_ws/src
cp -r ~/IKFoM/Samples/LINS---LiDAR-inertial-SLAM LINS---LiDAR-inertial-SLAM
cd ..
catkin_make -j1
source devel/setup.bash
```

**LIO-SAM Dataset:** Can be found [here](https://github.com/TixiaoShan/LIO-SAM), in which, Campus_dataset(large) and Campus_dataset(small) are tested in our paper.

Run
```
roslaunch lins run_port_exp.launch
rosbag play YOUR_DOWNLOADED.bag --clock
```

## 6.Acknowledgments
Thanks for C. Hertzberg,  R.  Wagner,  U.  Frese,  and  L.  Schroder.  Integratinggeneric   sensor   fusion   algorithms   with   sound   state   representationsthrough  encapsulation  of  manifolds.

