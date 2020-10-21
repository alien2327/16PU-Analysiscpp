# 16PU Analysis program by C/C++

このプログラムは16電極ビームモニターのデータ解析のためにC/C++で作成されている。

# Requirement

特に必要なパッケージはないが、g++コンパイラのバージョンが5.0以上、C++11の環境で一番スムーズに動く。

# Note

FFTやLinear fitting等、処理量が多い作業が含まれるため、従来のpythonスクリプトからこちらに移行したのである。  
また、データが入出力されるアドレスを直感的に確認することができるため、検証が終わった後にファームウェアに移植も用意になるかと思う。  
FPGAがC/C++の認識ができるので、直接入れるのもありかもしれない。

# Author

作成者　　:　李　耀漢（LEE YOHAN）  
所属　　　:　京都大学　理学研究科　高エネルギー研究室・J-PARCビームモニターチーム  
E-mail　　:　lee.yohan.83w@st.kyoto-u.ac.jp

# License

16PU system is under J-PARC intranet.

16PU system is Confidential.
