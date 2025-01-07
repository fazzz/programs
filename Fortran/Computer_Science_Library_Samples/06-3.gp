reset                       #設定のリセット
set view 0,0                #真上から（z軸に沿ってxy平面を）見るように視点を変更
set view equal xyz          #描画のスケールをxyz軸全て同じにする
set dgrid3d 100,100,1       #非格子状データから格子状データへの写像機能（x方向格子数,y方向格子数,1）
set contour                 #等値線をxy平面に描画
unset surface               #3次元プロットは非表示（等値線だけ描画する）
unset key                   #凡例を非表示
unset xtics                 #x軸の数値を非表示
unset ytics                 #y軸の数値を非表示
unset ztics                 #z軸の数値を非表示
splot '06-3b.txt' with line  #'06-3.txt'を描画（離散点じゃなくて線で）

# このファイルの使い方は「gnuplot> load '06-3.gp'」