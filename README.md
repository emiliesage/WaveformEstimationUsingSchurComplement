# WaveformEstimationUsingSchurComplement
## Introduction
The premise of the project is to create a method to estimate the phase, amplitude, offset, and frequency based on data from a noisy signal. The main method for estimating the signal is the least squares method. This is a robust method for estimating the amplitude, phase and offset of the signal when the frequency of the signal is known and the signal has introduced noise into the other variables of the signal. This method however is only useful if the frequency is known with great accuracy. when the frequency has noise and variation in itself, the estimations from the method become significantly less accurate. This is where the schur complement can help take into the account the error. The schur complement is a method that has multiple applications in the fields of statistics and linear algebra. In the statistical world, it takes multivariate gaussian random variables one of known data and one of varied data, and finds the expected value of the conditional data between the known and unknown. If we are able to find this expected value, it is possible to find a better estimation of the signal data based on a signal with noise in the frequency. this is the goal of the project. further discussion of the methods will be discussed in the methododlogy section.
## Materials
The materials used for this project include a physical and software component. The saoftware used in this project is Matlab, or in place of that, Octave, and Waveforms. The hardware used for this project is a Rhode Schwarz 4400 ociliscope, for the function generator, and an AD2 microcontroller with ADC (Analog to Digital Converter).
## Methodology
### Least Squares Method
Any signal with a single sinusoidal component can be described using an amplitude, 
frequency, phase, and offset, and can be written in the form

$$
r(t_n) = A \sin(\omega_0 t_n + \phi) + k
$$

When a signal is noisy, the linear least-squares sinusoidal method should be used to estimate the signal components $A$, $\phi$, and $k$. To find the parameters, the error 
must be defined. The error can be written as

$$
e(t_n) = s(n) - r(t_n)
$$

where $s(n)$ is the measured voltage. $r(t_n) = [E] \cdot p$, where the matrix $[E]$ can be defined as 

$$
[E] = \begin{bmatrix} \sin(\omega_0t_n) & \cos(\omega_0t_n) & 1 \end{bmatrix}
$$

and the $p$ vector written of the form

$$
\vec{p} = \begin{bmatrix} \alpha \\\\ \beta \\\\ k \end{bmatrix}
$$

where $\alpha = A \cos(\phi)$, $\beta = A \sin(\phi)$, and $k$ is the signal offset. Therefor to reconstruct a noisy signal in this model and find the amplitude, phase and offset of the estimated signal, the vector $p$ can be found by the formula:

$$
\vec{p} = ([E]^T[E])^{-1}[E]^T \cdot s(n)
$$

where $s(n)$ is the noisy signal data in the discrete time domain. This results in the $p$ vector shown above that discribes the signal parameters. Then the parameters can be calculated from this vector where $A = \sqrt{\alpha^2 + \beta^2}$, $\phi = \frac{\alpha}{\beta}$, and $k$. These can then be used to discribe an approximation to the signal from the discrete data.

## Results

![](https://github.com/emiliesage/WaveformEstimationUsingSchurComplement/figures/117khz.png)
