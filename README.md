# SPHCODE
Smoothed Particle Hydrodynamics (SPH)法のサンプルコードです。圧縮性流体専用です。

## 実装
### SPH方程式
#### Standard SPH (density-energy formulation; Springel & Hernquist 2002, Monaghan 2002)
最小作用の原理から導出したSPH法の方程式です。smoothing length の空間微分 (grad-h term)を考慮した方程式系になっています。

#### Density Independent SPH (pressure-energy formulation; Saitoh & Makino 2013; Hopkins 2013)
状態方程式からSPH粒子の体積を求めることによって、密度に陽に依存しない方程式にします。接触不連続面を正しく扱えるようになります。

#### Godunov SPH (Inutsuka 2002; Cha & Whitworth 2003; Murante et al. 2011)
Riemann solverを使って粒子間相互作用を計算することで、人工粘性のような計算を安定化させる数値拡散が自動的に入るようになります。

### カーネル関数
#### Cubic spline (e.g. Monaghan 1992)
昔からよく使われているオーソドックスなカーネル関数です。

#### Wendland C4 (Wendland 1995)
Cubic splineカーネルより高次のカーネル関数です。影響半径内の粒子数が大きいときに粒子同士がくっついてしまう不安定性 (pairing instability; Dehnen & Aly 2012) を防ぐことができます。

### 人工粘性
#### signal velocity formulation (Monaghan 1997)
Riemann solverによる数値拡散から類推して導出された人工粘性です。

#### Balsara switch (Balsara 1995)
速度の回転が発散より大きい領域で人工粘性係数を小さくすることで、シアー領域に余分な粘性が効かないようにします。

#### 時間依存人工粘性係数 (Rosswog et al. 2000)
SPH粒子それぞれの人工粘性係数を時間変化させます。圧縮領域では人工粘性係数を大きく、それ以外では小さくするようにします。

### その他
#### 人工熱伝導 (Price 2008; Wadsley et al. 2008)
#### 自己重力 (Hernquist & Katz 1989)
#### Tree (Hernquist & Katz 1989)


## 参考文献
* Balsara, D. S. (1995). von Neumann stability analysis of smoothed particle hydrodynamics--suggestions for optimal algorithms. Journal of Computational Physics, 121(2), 357–372. https://doi.org/10.1016/S0021-9991(95)90221-X
* Cha, S. H., & Whitworth, A. P. (2003). Implementations and tests of Godunov-type particle hydrodynamics. Monthly Notices of the Royal Astronomical Society, 340(1), 73–90. https://doi.org/10.1046/j.1365-8711.2003.06266.x
* Dehnen, W., & Aly, H. (2012). Improving convergence in smoothed particle hydrodynamics simulations without pairing instability. Monthly Notices of the Royal Astronomical Society, 425(2), 1068–1082. https://doi.org/10.1111/j.1365-2966.2012.21439.x
* Hernquist, L., & Katz, N. (1989). TREESPH - A unification of SPH with the hierarchical tree method. The Astrophysical Journal Supplement Series, 70, 419. https://doi.org/10.1086/191344
* Hopkins, P. F. (2013). A general class of Lagrangian smoothed particle hydrodynamics methods and implications for fluid mixing problems. Monthly Notices of the Royal Astronomical Society, 428(4), 2840–2856. https://doi.org/10.1093/mnras/sts210
* Inutsuka, S. (2002). Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver. Journal of Computational Physics, 179, 238. https://doi.org/10.1006/jcph.2002.7053
* Murante, G., Borgani, S., Brunino, R., & Cha, S.-H. (2011). Hydrodynamic simulations with the Godunov smoothed particle hydrodynamics. Monthly Notices of the Royal Astronomical Society, 417(1), 136–153. https://doi.org/10.1111/j.1365-2966.2011.19021.x
* Monaghan, J. J. (1992). Smoothed Particle Hydrodynamics. Annual Review of Astronomy and Astrophysics, 30(1), 543–574. https://doi.org/10.1146/annurev.aa.30.090192.002551
* Monaghan, J. J. (1997). SPH and Riemann Solvers. Journal of Computational Physics, 136(2), 298–307. https://doi.org/10.1006/jcph.1997.5732
* Monaghan, J. J. (2002). SPH compressible turbulence. Monthly Notices of the Royal Astronomical Society, 335(3), 843–852. https://doi.org/10.1046/j.1365-8711.2002.05678.x
* Price, D. J. (2008). Modelling discontinuities and Kelvin–Helmholtz instabilities in SPH. Journal of Computational Physics, 227(24), 10040–10057. https://doi.org/10.1016/j.jcp.2008.08.011
* Rosswog, S., Davies, M. B., Thielemann, F.-K., & Piran, T. (2000). Merging neutron stars: asymmetric systems. Astronomy and Astrophysics, 360, 171–184. Retrieved from http://adsabs.harvard.edu/abs/2000A&A...360..171R%5Cnpapers2://publication/uuid/39C9D6F4-C091-499D-8F66-867A98C4DD32
* Saitoh, T. R., & Makino, J. (2013). A DENSITY-INDEPENDENT FORMULATION OF SMOOTHED PARTICLE HYDRODYNAMICS. The Astrophysical Journal, 768(1), 44. https://doi.org/10.1088/0004-637X/768/1/44
* Springel, V., & Hernquist, L. (2002). Cosmological smoothed particle hydrodynamics simulations: the entropy equation. Monthly Notices of the Royal Astronomical Society, 333(3), 649–664. https://doi.org/10.1046/j.1365-8711.2002.05445.x
* Wadsley, J. W., Veeravalli, G., & Couchman, H. M. P. (2008). On the treatment of entropy mixing in numerical cosmology. Monthly Notices of the Royal Astronomical Society, 387(1), 427–438. https://doi.org/10.1111/j.1365-2966.2008.13260.x
* Wendland, H. (1995). Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree. Advances in Computational Mathematics, 4(1), 389–396. https://doi.org/10.1007/BF02123482
