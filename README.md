# WaveformEstimationUsingSchurComplement
## Introduction
The premise of the project is to create a method to estimate the phase, amplitude, offset, and frequency based on data from a noisy signal. The main method for estimating the signal is the least squares method. This is a robust method for estimating the amplitude, phase and offset of the signal when the frequency of the signal is known and the signal has introduced noise into the other variables of the signal. This method however is only useful if the frequency is known with great accuracy. when the frequency has noise and variation in itself, the estimations from the method become significantly less accurate. This is where the schur complement can help take into the account the error. The schur complement is a method that has multiple applications in the fields of statistics and linear algebra. In the statistical world, it takes multivariate gaussian random variables one of known data and one of varied data, and finds the expected value of the conditional data between the known and unknown. If we are able to find this expected value, it is possible to find a better estimation of the signal data based on a signal with noise in the frequency. this is the goal of the project. further discussion of the methods will be discussed in the methododlogy section.
## Materials
The materials used for this project include a physical and software component. The saoftware used in this project is Matlab, or in place of that, Octave, and Waveforms. The hardware used for this project is a Rhode Schwarz 4400 ociliscope, for the function generator, and an AD2 microcontroller with ADC (Analog to Digital Converter).
## Methodology
### Least Squares Method
Any signal with a single sinusoidal component can be described using an amplitude, 
frequency, phase, and offset, and can be written in the form [2]

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

where $s(n)$ is the noisy signal data in the discrete time domain. This results in the $p$ vector shown above that discribes the signal parameters. Then the parameters can be calculated from this vector where $A = \sqrt{\alpha^2 + \beta^2}$, $\phi = \frac{\alpha}{\beta}$, and $k$. These can then be used to discribe an approximation to the signal from the discrete data. This method is very accurate and robust when the frequency of the signal is known and has no variation but when there is variation in the frequency then the estimation of amplitude and phase drops off significantly. This is why the schur complement is useful.

### Schur Complement
The Schur complement is a linear algebra form that is used in the areas of statistics. In statistical application, the goal of this model is to find the conditional expected value of a specific set of multivariant data in relation to another dataset, taking into account the mean and covariance of the set.

The first step to seting up this model is to assume that the datasets are able to be modeled as a gaussian random variable. This means that the multivariant system can be expressed as the equation:

$$
	\cal{N}_k(\mu,\Sigma) \frac{1}{\sqrt{(2 \pi)^k|\Sigma|}}e^{(-\frac{(x - \mu)^T \Sigma^{-1}(x - \mu)}{2})}
$$

Where $\mu$ is the mean of the multivariant system, $k$ is the number of systems, and $\Sigma$ is the covariance matrix showing the relationship between each variable in the system. so long as the system can be estimated in this way the schur complement method can be used. 

The first step of this process was to run a Monte-Carlo simulation to train the vectors to be used in the Schur complement. A signal was generated for $5 \cdot 10^3$ simulations to tran the datset. For one random variable called $X$ the amplitude, phase, offset and frequency were calculated with noise and direct values were given for each. This is to be the reference set. Then noise was added to the signal and the least squares approximation was evaluated for this signal. This was stored in the variable $Y$. This is the varied set used as a random variable in the calculation. For frequncy estimation in this reguard, The fast forier transform of the signal was also used in the approximation for use in live data.

Then the two sets were coupled into a random variable 

$$
z = \begin{bmatrix} x \\\\ y \end{bmatrix} \textasciitilde \cal{N} \begin{pmatrix} \begin{bmatrix}\mu_x \\\\ \mu_y \end{bmatrix} & \begin{bmatrix} \Sigma_{xx} & \Sigma_{xy} \\\\ \Sigma_{yx} & \Sigma_{yy} \end{bmatrix} \end{pmatrix}
$$

From this random variable we can use schur complement to find the conditional covariance and also the conditional expected value. These are of the form [1]:

$$
\mu_{x|y} = \mu_x + \Sigma_{xy} \Sigma_{yy}^{-1} (y' - \mu_y), ~ \Sigma_{x|y} = \Sigma_{xx} - \Sigma_{xy} \Sigma_{yy}^{-1} \Sigma_{yx}
$$

To derive the formula the conditional covariance $\Sigma_{x|y}$ naturally arives through taking the schur complement of the joint covariance matrix. The following derivation is the completing the square method to derive the conditional expected value $\mu_{x|y}$.

The conditional joint density can be written of the form:

$$
	p(x,y) ~ \alpha ~ exp(-\frac{1}{2} \begin{pmatrix}z \\\\ d \end{pmatrix}^T \Sigma^{-1} \begin{pmatrix}z \\\\ d \end{pmatrix})
$$

where the variables are:

$$
	z = x - \mu_x, ~ d = y - \mu_y, ~ \Sigma^{-1} = \begin{pmatrix} M & N \\\\ N^T & P \end{pmatrix}
$$

The conditional probability after the multiplication is completed becomes:

$$
	p(x|y) ~ \alpha ~ exp(-\frac{1}{2}(z^TMz + 2z^TNd))
$$

By completinf the square we can reevaluate the equation to look like this

$$
	p(x|y) ~ \alpha ~ exp(\frac{1}{2}(z + M^{-1}Nd)^T M(z + M^{-1}Nd)) \times (constant in y)
$$

## Results

![118khz signal](https://github.com/emiliesage/WaveformEstimationUsingSchurComplement/blob/main/figures/118khz.png)

## Conclusion

## References

[1] Santos, Talles, DEVELOPMENT AND USE OF ANATOMICAL AND
PHYSIOLOGICAL PRIOR INFORMATION TO
ESTIMATE ELECTRICAL IMPEDANCE
TOMOGRAPHY IMAGE, 2019

[2] Lombello, Christiane Bertachini, and Patricia Aparecida da Ana. Current Trends in Biomedical Engineering. Springer, 2023.
