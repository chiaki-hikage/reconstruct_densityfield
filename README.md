# 内容
本コードは、物質(ダークマター)の初期位置からの変位を戻すことで、
重力成長した物質密度場から線形成長したゆらぎを回復するコードです。

宇宙の物質は重力によって集まって成長し、現在観測される大規模な構造ができます。しかし重力成長は非線形な効果であるため、
異なるスケールのゆらぎが互いに相関をもつようになり、精密な理論モデルを構築すること難しくなります。
またゆらぎの統計がガウス統計から大きくずれるため、宇宙論的情報を十分に引き出すのにより高次の統計量を使う必要が出てきます。

ラグランジュ的な記述では物質素片の運動は初期位置からの変位ベクトルで記述されます。
初期位置からの変位のうち、大スケールの速度場(バルクモーション)に由来する成分は、
重力成長後の密度場から近似的に求めることができ、その影響を補正することで初期の線形に近いゆらぎに戻すことができます。
この操作によって、重力成長の非線形性の影響を部分的に取り除くことができ、解析が容易になります。

この方法は、もともとは宇宙大規模構造に含まれるバリオン音響振動のシグナルを改善するために
用いられたものですが(Eistenstein et al. 2007)、非線形な影響が補正される効果で、
摂動論の適用範囲が伸びたり、高次統計量に流出した情報がパワースペクトルのもつ宇宙論的情報が増えるなど
さまざまな恩恵があることが分かりました。

# References
- Perturbation Theory for BAO reconstructed fields: one-loop results in real-space matter density field  
Chiaki Hikage, Kazuya Koyama, Alan Heavens  
Phys. Rev. D, Vol. 96 (2017) id.043513

Perturbation theory for the redshift-space matter power spectra after reconstruction  
Chiaki Hikage, Kazuya Koyama, Ryuichi Takahashi  
Phys. Rev. D, Vol. 101, Issue 4 (2020), id.043510
