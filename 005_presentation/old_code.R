## -- old code mod 3 presentation 
<!--
        
.footnote[Gauch & Whittaker (1972), Zuur et al. (2007)]
# The idea is old 

The Gaussian response model: 
        $$Y_{is} = c_s\ exp\bigg(- \frac{(x_{ip} - u_{sp})^2}{2t_{sp}^2}\bigg) $$
        
        .blue[c]: maximal abundance, .blue[u]: location of optima, .blue[t]: tolerance  

--
        
        As GLM:
        \begin{align}
y_{is} &= exp \bigg(ln(c_s) - \frac{u_{sp}^2}{2t_{sp}^2} + \frac{u_{sp}}{t_{sp}^2}x_{ip} - \frac{1}{2t_{sp}^2}x^2_{ip}\bigg)\\
&= exp(b_{1s} + b_{sp}\ x_{ip} + b_{sp}\ x_{ip}^2)
\end{align}

with  

$$t_s = \frac{1}{\sqrt{-2b_{sp}}}; u_{sp} = \frac{-b_{sp}}{2b_{sp}}; c_s = exp \bigg(b_{sp} - \frac{b_{sp}^2}{4b_{sp}}\bigg)$$ 
        
        ---
        # Restricted Gaussian Regression 
        
        `r icon::fa("cog")``r icon::fa("cog", color = "white")` Uses $(1+2P)S$ parameters.  

--
        
        `r icon::fa("chevron-circle-up")``r icon::fa("cog", color = "white")`For 10 species and 5 variables 110 parameters.   

--
        
        Gradients $\mathbf{\Psi}$ as linear combinations of measured variables $x$:
        $$\psi_r = \Sigma_{p =  1}^P \alpha_{rp}\ x_p$$
        --
        
        Plug gradients back into Gaussian regression
$$Y_{is} = c_s\ exp\bigg(\frac{(\psi_{ir} - u_{sr})^2}{2t_{sr}^2}\bigg)$$
        
        If R $\ll$ P, this reduces the number of parameters. 
-->