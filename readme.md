## Inverted pendulum for simulation

We use an inverted pendulum for simulation, the physical parameters are shown in the following figure:

<img src="/figs/inv_pen.png" width="480">

Based on the knowledge of theoretical mechanics, one obtains the inverted pendulum system dynamic:

<a href="https://www.codecogs.com/eqnedit.php?latex=(M&plus;m)\ddot{x}&plus;b&space;\dot{x}&space;&plus;ml&space;\ddot{\theta}&space;\cos&space;\theta&space;-&space;ml&space;\dot{\theta}^2\sin\theta=F" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(M&plus;m)\ddot{x}&plus;b&space;\dot{x}&space;&plus;ml&space;\ddot{\theta}&space;\cos&space;\theta&space;-&space;ml&space;\dot{\theta}^2\sin\theta=F" title="(M+m)\ddot{x}+b \dot{x} +ml \ddot{\theta} \cos \theta - ml \dot{\theta}^2\sin\theta=F" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=ml^2&space;\ddot{\theta}&space;&plus;mgl\sin\theta&space;=ml\ddot{x}\cos\theta." target="_blank"><img src="https://latex.codecogs.com/gif.latex?ml^2&space;\ddot{\theta}&space;&plus;mgl\sin\theta&space;=ml\ddot{x}\cos\theta." title="ml^2 \ddot{\theta} +mgl\sin\theta =ml\ddot{x}\cos\theta." /></a>


Rewrite it as a non-linear state space model :

<a href="https://www.codecogs.com/eqnedit.php?latex=\left[\begin{array}{cccc}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;m&plus;M&space;&&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;&&space;0&space;&&space;m&space;l^{2}&space;\end{array}\right]&space;\frac{d}{d&space;t}\left[\begin{array}{c}&space;x&space;\\&space;\dot{x}&space;\\&space;\theta&space;\\&space;\dot{\theta}&space;\end{array}\right]=\left[\begin{array}{c}&space;\dot{x}&space;\\&space;m&space;l&space;\dot{\theta}^{2}&space;\sin&space;\theta-b&space;\dot{x}&space;\\&space;\dot{\theta}&space;\\&space;m&space;g&space;l&space;\sin&space;\theta&space;\end{array}\right]&plus;\left[\begin{array}{c}&space;0&space;\\&space;1&space;\\&space;0&space;\\&space;0&space;\end{array}\right]&space;F." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[\begin{array}{cccc}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;m&plus;M&space;&&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;&&space;0&space;&&space;m&space;l^{2}&space;\end{array}\right]&space;\frac{d}{d&space;t}\left[\begin{array}{c}&space;x&space;\\&space;\dot{x}&space;\\&space;\theta&space;\\&space;\dot{\theta}&space;\end{array}\right]=\left[\begin{array}{c}&space;\dot{x}&space;\\&space;m&space;l&space;\dot{\theta}^{2}&space;\sin&space;\theta-b&space;\dot{x}&space;\\&space;\dot{\theta}&space;\\&space;m&space;g&space;l&space;\sin&space;\theta&space;\end{array}\right]&plus;\left[\begin{array}{c}&space;0&space;\\&space;1&space;\\&space;0&space;\\&space;0&space;\end{array}\right]&space;F." title="\left[\begin{array}{cccc} 1 & 0 & 0 & 0 \\ 0 & m+M & 0 & m l \cos \theta \\ 0 & 0 & 1 & 0 \\ 0 & m l \cos \theta & 0 & m l^{2} \end{array}\right] \frac{d}{d t}\left[\begin{array}{c} x \\ \dot{x} \\ \theta \\ \dot{\theta} \end{array}\right]=\left[\begin{array}{c} \dot{x} \\ m l \dot{\theta}^{2} \sin \theta-b \dot{x} \\ \dot{\theta} \\ m g l \sin \theta \end{array}\right]+\left[\begin{array}{c} 0 \\ 1 \\ 0 \\ 0 \end{array}\right] F." /></a>,
where $b$ is the coefficient of friction for cart movement and is set as $0$ in the simulation.

Reformulate the system as <a href="https://www.codecogs.com/eqnedit.php?latex=\dot{x}=f(x)&plus;g(x)u" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dot{x}=f(x)&plus;g(x)u" title="\dot{x}=f(x)+g(x)u" /></a> and linearize system at <a href="https://www.codecogs.com/eqnedit.php?latex=\theta=0,\dot{\theta}=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta=0,\dot{\theta}=0" title="\theta=0,\dot{\theta}=0" /></a>.
Then one obatins a continuous time LTI system: $\dot{x}(t)=A_c x(t) +B_c u(t)$.
Sample the continuous time system with period $T_s$, then the discrete system is
$x(k+1)=Ax(k)+Bu(k)$ where

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." title="A=\exp(T\cdot A_c),\ B=\left(\int_{0}^{T} \exp(\tau A_c) d\tau \right)B_c ." /></a>

The system output equation is

<a href="https://www.codecogs.com/eqnedit.php?latex=Y(k)=Cx(k)&plus;v(k)&plus;a(k),\text{&space;where&space;}&space;C=\begin{bmatrix}1&0&0&0&space;\\\&space;1&0&0&0&space;\\\&space;1&0&0&0&space;\\\&space;0&0&1&0&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y(k)=Cx(k)&plus;v(k)&plus;a(k),\text{&space;where&space;}&space;C=\begin{bmatrix}1&0&0&0&space;\\\&space;1&0&0&0&space;\\\&space;1&0&0&0&space;\\\&space;0&0&1&0&space;\end{bmatrix}" title="Y(k)=Cx(k)+v(k)+a(k),\text{ where } C=\begin{bmatrix}1&0&0&0 \\\ 1&0&0&0 \\\ 1&0&0&0 \\\ 0&0&1&0 \end{bmatrix}" /></a>

where $v(k)$ is the output noise and $a(k)$ is the injected attack.
There are 3 sensors monitoring the cart position and 1 sensor monitoring the pendulum angle.

The system is simulated with sampling time $T_s=0.1 s$.

The system noise process has covariance <a href="https://www.codecogs.com/eqnedit.php?latex=Q=T_s^2&space;\cdot&space;{\rm&space;diag}[0.1\&space;0.1\&space;0.01\&space;0.01]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q=T_s^2&space;\cdot&space;{\rm&space;diag}[0.1\&space;0.1\&space;0.01\&space;0.01]" title="Q=T_s^2 \cdot {\rm diag}[0.1\ 0.1\ 0.01\ 0.01]" /></a>
and the measurement noise is <a href="https://www.codecogs.com/eqnedit.php?latex=R=T_s^2&space;\cdot&space;{\rm&space;diag}[0.1\&space;0.1\&space;0.1\&space;0.1]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R=T_s^2&space;\cdot&space;{\rm&space;diag}[0.1\&space;0.1\&space;0.1\&space;0.1]" title="R=T_s^2 \cdot {\rm diag}[0.1\ 0.1\ 0.1\ 0.1]" /></a>,
i.e., the noise is scaled according to sampling time $T_s$. The initial state is set as $x_0=[0,1,0,1]'$, and the estiamtor is assumed to know the initial state.

### Remark on state transformation

In order to better analyze the observability structure, the matrix $A$ is transformed into a Jordan canonical form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\bar{A}=U^{-1}AU=\begin{bmatrix}&space;1&space;&&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0.642&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1.557&space;\end{bmatrix},\text{&space;where&space;}&space;U=\begin{bmatrix}&space;0.050&space;&&space;0.505&space;&&space;-0.253&space;&&space;-0.252\\&space;0&space;&&space;0.500&space;&&space;1.108&space;&&space;-1.108\\&space;0&space;&&space;0&space;&&space;0.500&space;&&space;0.500&space;\\&space;0&space;&&space;0&space;&&space;-2.217&space;&&space;2.217&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bar{A}=U^{-1}AU=\begin{bmatrix}&space;1&space;&&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0.642&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1.557&space;\end{bmatrix},\text{&space;where&space;}&space;U=\begin{bmatrix}&space;0.050&space;&&space;0.505&space;&&space;-0.253&space;&&space;-0.252\\&space;0&space;&&space;0.500&space;&&space;1.108&space;&&space;-1.108\\&space;0&space;&&space;0&space;&&space;0.500&space;&&space;0.500&space;\\&space;0&space;&&space;0&space;&&space;-2.217&space;&&space;2.217&space;\end{bmatrix}" title="\bar{A}=U^{-1}AU=\begin{bmatrix} 1 & 1 & 0 & 0\\ 0 & 1 & 0 & 0\\ 0 & 0 & 0.642 & 0 \\ 0 & 0 & 0 & 1.557 \end{bmatrix},\text{ where } U=\begin{bmatrix} 0.050 & 0.505 & -0.253 & -0.252\\ 0 & 0.500 & 1.108 & -1.108\\ 0 & 0 & 0.500 & 0.500 \\ 0 & 0 & -2.217 & 2.217 \end{bmatrix}" /></a>

And the simulation is performed on the following system where $\bar{x}=U^{-1}x$:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;\bar{x}(k)&=&space;\bar{A}&space;\bar{x}(k)&space;&plus;&space;U^{-1}B\bar{x}(k)&plus;U^{-1}w(k),&space;\\&space;y(k)&=CU\bar{x}(k)&plus;v(k)&plus;a(k)&space;,&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;\bar{x}(k)&=&space;\bar{A}&space;\bar{x}(k)&space;&plus;&space;U^{-1}B\bar{x}(k)&plus;U^{-1}w(k),&space;\\&space;y(k)&=CU\bar{x}(k)&plus;v(k)&plus;a(k)&space;.&space;\end{align*}" title="" /></a>

The estimation of $x(k)$ can be obtained by performing transformation $U$ on the estimation of $\hat{x}(k)$.
