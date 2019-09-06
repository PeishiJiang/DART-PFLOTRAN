missingvals = -888888.0000;
nx =        1;
ny =        5;
nz =        7;
nens =        1;
interptest = [ ... 
  0.11957206792040614     
  0.10447324792485997     
  0.10533243053066269     
  0.10173719412599454     
  0.10201132211106431     
  0.10081335999019037     
  0.10092189663355389     
  0.11927507566385899     
  0.10447324792485999     
  0.10526554103808819     
  0.10173719412599454     
  0.10198998056888466     
  0.10081335999019034     
  0.10091344678721006     
  0.11957206792040614     
  0.10447324792485997     
  0.10533243053066269     
  0.10173719412599454     
  0.10201132211106431     
  0.10081335999019037     
  0.10092189663355389     
  0.11927507566385899     
  0.10447324792485999     
  0.10526554103808819     
  0.10173719412599454     
  0.10198998056888466     
  0.10081335999019034     
  0.10091344678721006     
  0.11957206792040614     
  0.10447324792485997     
  0.10533243053066269     
  0.10173719412599454     
  0.10201132211106431     
  0.10081335999019037     
  0.10092189663355389     
];
datmat = reshape(interptest,nz,ny,nx,nens);
datmat = permute(datmat,[4,1,2,3]);
datmat(datmat == missingvals) = NaN;
