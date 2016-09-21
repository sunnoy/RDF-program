pro THIIS ;,xrdfile,bgfile,u,rho0,c1,c2,c3,c4,t,R,S,lambda,alpha
  ;data input begin
  ;get PDF from XRD data,(for example:xrdfile='C:\Users\IDLWorkspace82\get\xrd.txt')
  ;input the information and parameters
  xrdfile=
  bgfile=
  u=0.32 ;refer to the DAT file: 'Um - atomic weight.dat'; rho is g/cm^3
  rho0=0.35 ;the value of
  ;Scattering factors: 1-C,2-H,3-N,4-0;
  c1=0.25 & c2=0.35 & c3=0.25 & c4=0.36; atomic concentration;
  t=0.214;mm
  R=0.351;mm
  S=0.213;rad
  lambda=0.701; wavelength - Mo=0.701; Cu=1.54
  alpha=0.014
  ;data input end
  
  
  ;Cromer parameters: for atomic factor, were not affected by Mo or Cu target
  cromer1=[2.31,20.8439,1.02,10.2075,1.5886,0.5687,0.865,51.6512,0.2156]
  cromer2=[0.489918,20.6593,0.262003,7.74039,0.196767,49.5519,0.049879,2.20159, 0.001305]
  cromer3=[12.2126,0.0057,3.1322,9.8933,2.0125,28.9975,1.1663,0.5826,-11.529]
  cromer4=[3.0485,13.2771,2.2868,5.7011,1.5463,0.3239,0.867,32.9089,0.2508]
  ;Incoherently Scattered X-Ray Intensities:
  ; method 1: balyuzi
  ;; the sequence of balyuzi=a1, b1,a2, b2,a3,b3, a4, b4, a5,b5; were not affected by Mo or Cu target
  balyuzi1=[0.756799996,82.2385025,2.55110002,31.7282009,0.7051,11.9470997,1.4605,1.46370006,0.526300013,0.514999986]
  balyuzi2=[0.262300014,32.3717003,0.509400010,14.7083998,0.203400001,6.68839979,0.02490000054,2.48429990,0.0,0.0]
  balyuzi3=[0.907000005,64.1555023,2.89720011,20.8507004,1.16589999,7.75759983,1.55260003,1.03349996,0.476900011,0.351599991]
  balyuzi4=[0.8847,52.0063019,3.21889997,16.4487,1.79900002,6.59579992,1.55379999,0.814300001,0.54339999,0.281500012]
  ;read data from files.
  xrdRead=read_ascii(xrdfile)
  ;read backgroudfile.
  bgRead=read_ascii(bgfile)
  im=xrdRead.(0) & ib=bgRead.(0)
  theta2=im(0,*)*!pi/180.0 & theta=theta2/2.0
  ar=0.5+(0.5-t*cos(theta)/(R*S))*exp(-2.0*u*t/sin(theta))
  ;at=(1.0-t*sin(theta)/(R*S))*exp(-2.0*u*t/cos(theta))
  ia=ib & ia(1,*)=ib(1,*)*ar
  Atheta=(1.0-exp(-2*u*t/sin(theta)))/(2.0*u)
  Ptheta=(1.0+0.9563*(cos(theta2))^2)/2.0
  ic=im & ic(1,*)=(im(1,*)-ia(1,*))/(Atheta*Ptheta)
  ;2Theta to Q.
  iq=ic
  iq(0,*)=4.0*!pi*sin(theta)/lambda & iq(1,*)=smooth(iq(1,*),13,/edge_truncate)
  Qmax=max(iq(0,*)) & Qmin=min(iq(0,*)) & deltaQ=0.02
  Q=[0.0] & for Qi=Qmin,Qmax,deltaQ do Q=[Q,Qi] & Q=Q(1:*)
  ipy=interpol(iq(1,*),iq(0,*),Q,/lsquadratic) & ip=transpose([[Q],[ipy]])
  ;Incoherently Scattered X-Ray Intensities: Compton scattering;
  Qp=(Q/(4.0*!pi))^2
  ;Balyuzi Analytic Approximation
  in1=0.0 & in2=0.0 & in3=0.0 & in4=0.0
  for i=1,10,2 do begin
    ai=balyuzi1(i-1) & bi=balyuzi1(i)
    in1 =in1 +ai*exp(-bi*Qp)
  endfor
  for i=1,10,2 do begin
    ai=balyuzi2(i-1) & bi=balyuzi2(i)
    in2=in2+ai*exp(-bi*Qp)
  endfor
  for i=1,10,2 do begin
    ai=balyuzi3(i-1) & bi=balyuzi3(i)
    in3=in3+ai*exp(-bi*Qp)
  endfor
  for i=1,10,2 do begin
    ai=balyuzi4(i-1) & bi=balyuzi4(i)
    in4=in4+ai*exp(-bi*Qp)
  endfor
  in=c1*in1+c2*in2+c3*in3+c4*in4
  ;Cromer function, atomic scattering factor,
  ;fq=c+sum(ai*exp(-bi*((Q^2)/(16*!pi^2)), i=from 1 to 4
  fq1=cromer1(8) & fq2=cromer2(8) & fq3=cromer3(8) & fq4=cromer4(8)
  for i=1,8,2 do begin
    ai=cromer1(i-1) & bi=cromer1(i)
    fq1=fq1+ai*exp(-bi*(Q^2)/(16*!pi^2))
  endfor
  for i=1,8,2 do begin
    ai=cromer2(i-1) & bi=cromer2(i)
    fq2=fq2+ai*exp(-bi*(Q^2)/(16*!pi^2))
  endfor
  for i=1,8,2 do begin
    ai=cromer3(i-1) & bi=cromer3(i)
    fq3=fq3+ai*exp(-bi*(Q^2)/(16*!pi^2))
  endfor
  for i=1,8,2 do begin
    ai=cromer4(i-1) & bi=cromer4(i)
    fq4=fq4+ai*exp(-bi*(Q^2)/(16*!pi^2))
  endfor
  f2a=c1*fq1^2+c2*fq2^2+c3*fq3^2+c4*fq4^2
  fa2=(c1*fq1+c2*fq2+c3*fq3+c4*fq4)^2
  ;normalization factor: beta.
  ;RDF method for calculating beta factor.
  numQ=n_elements(Q) & beta1=0.0 & beta2=0.0
  for i=1,numQ do begin
    Qi=ip(0,i-1) & ini=in(i-1) & f2ai=f2a(i-1) & fa2i=fa2(i-1)
    beta1=beta1+(f2ai+ini)/fa2i*Qi^2*exp(-Qi^2*alpha^2)*deltaQ
    beta2=beta2+ip(1,i-1)/fa2i*Qi^2*exp(-Qi^2*alpha^2)*deltaQ
  endfor
  normBeta=(beta1-2.0*!pi^2*rho0)/beta2
  ;high angle method for calculating factor.
  ;numQ=n_elements(Q) & beta1=0.0 & beta2=0.0
  ;for i=1,numQ do begin
  ;Qi=ip(0,i-1) & ini=in(i-1) & f2ai=f2a(i-1) & fa2i=fa2(i-1)
  ;beta1=beta1 +(f2ai+ini)*deltaQ
  ;beta2=beta2+ip(1,i-1)/fa2i*deltaQ
  ;endfor
  ;normBeta=beta1/beta2
  ;calculate I(Q).
  is=ip & is(1,*)=(normBeta*ip(1,*)-f2a+fa2-in)/fa2
  ;FFT for G(r)
  Gr=[0.0,0.0] & rdf=[0.0,0.0] & averageLine=[0.0,0.0] & zeroLine=[0.0,0.0]
  for r=0.75,20.0,0.02 do begin
    Gri=0.0
    for i=1,numQ do begin
      Qi=is(0,i-1) & Sqi=Qi*(is(1,i-1)-1)
      Gri=Gri+2.0/!pi*Sqi*exp(-Qi^2*alpha^2)*sin(Qi*r)*deltaQ
    endfor
    Gr=[[Gr],[r,Gri]]
    rdf=[[rdf],[r,4*!pi*r^2*rho0+r*Gri]]
    if r lt 10.0 then averageLine=[[averageLine],[r,-4*!pi*r*rho0]]
    zeroLine=[[zeroLine],[r,0.0]]
  endfor
  GrData=Gr(*,1:*)
  pdfData=GrData & pdfData(1,*)=1.0+GrData(1,*)/(4.0*!pi*GrData(0,*)*rho0)
  rdfData=rdf(*, 1: *)
  ;plot G(r)
  plot,GrData(0,*),GrData(1,*)
  plots,zeroLine(0,1 :*),zeroLine(1 ,1 :*),color=255
  plots,averageLine(0,1 :*),averageLine(1,1:*),color=255
  ;save files.
  sn=strlen(xrdfile)
  IqPaths=strmid(xrdfile,0,sn-4)+'-iq.dat'
  openw,lun,IqPaths,/get_lun & printf,lun,is,format=' (f10.2,f15.3) ' & free_lun,lun
  SqPaths=strmid(xrdfile,0,sn-4)+ '-sq.dat' & SqData=is &
  SqData(1,*)=Q*(is(1,*)-1)
  openw,lun,SqPaths,/get_lun & printf,lun,SqData,format=' (f10.2,f15.3) '&
  free_lun,lun
  GrPaths=strmid(xrdfile,0,sn-4)+ '-rrdf.dat'
  openw,lun,GrPaths,/get_lun & printf,lun,GrData,format= ' (f10.3,f15.3) '&
  free_lun,lun
  pdfPaths=strmid(xrdfile,0,sn-4)+ '-pdf.dat '
  openw,lun,pdfPaths,/get_lun & printf,lun,pdfData,format='(f10.3,f15.3)' &
  free_lun,lun
  rdfPaths=strmid(xrdfile,0,sn-4)+ '-rdf.dat '
  openw,lun,rdfPaths,/get_lun & printf,lun,rdfData,format= ' (f10.3,f15.3) ' &
  free_lun,lun
  print ,'Calculation done and PDF data files has been calculated '
  
END

