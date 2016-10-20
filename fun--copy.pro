PRO TLB_EVENT, ev
  WIDGET_CONTROL, ev.TOP,get_uvalue= gz
  
  WIDGET_CONTROL, ev.id, GET_UVALUE=eventval
  case eventval of
  
    'xrd': BEGIN
      md = DIALOG_PICKFILE(filter ='*.txt')
      widget_control, gz.xrd,  SET_VALUE=md
      COMMON sss, ss, ssd ;设定全局变量
      ss = md   ;转换为全局变量
    END
    
    'xrdd': BEGIN
      mdd = DIALOG_PICKFILE(filter ='*.txt')
      widget_control, gz.xrdd,  SET_VALUE=mdd
      ssd = mdd
    END
    
    ;清空输入数据
    'reset': BEGIN
      widget_control, gz.u,  SET_VALUE=''
      widget_control, gz.p,  SET_VALUE=''
      widget_control, gz.c,  SET_VALUE=''
      widget_control, gz.h,  SET_VALUE=''
      widget_control, gz.n,  SET_VALUE=''
      widget_control, gz.o,  SET_VALUE=''
      widget_control, gz.t,  SET_VALUE=''
      widget_control, gz.r,  SET_VALUE=''
      widget_control, gz.s,  SET_VALUE=''
      widget_control, gz.l,  SET_VALUE=''
      widget_control, gz.a,  SET_VALUE=''
      widget_control, gz.xrd,  SET_VALUE=''
      widget_control, gz.xrdd,  SET_VALUE=''
    END
    
    ;启动起算
    'run': BEGIN
      aa = ss
      b = ssd
      widget_control, gz.u,  get_value=us
      u=FLOAT(us)
      widget_control, gz.p,  get_value=ps
      p=fLOAT(ps)
      widget_control, gz.c,  get_value=cs
      c=fLOAT(cs)
      widget_control, gz.h,  get_value=hs
      h=fLOAT(hs)
      widget_control, gz.n,  get_value=ns
      n=fLOAT(ns)
      widget_control, gz.o,  get_value=os
      o=fLOAT(os)
      widget_control, gz.t,  get_value=ts
      t=fLOAT(ts)
      widget_control, gz.r,  get_value=rs
      r=fLOAT(rs)
      widget_control, gz.s,  get_value=ss
      s=fLOAT(ss)
      widget_control, gz.l,  get_value=ls
      l=fLOAT(ls)
      widget_control, gz.a,  get_value=as
      a=fLOAT(as)
      
      pri, aa, b, u, p, c, h, n, o, t, r, s, l, a ;输入数据传入计算过程
      
    END
    
    
  endcase
  
END

pro pri, aa, b, u, p, c, h, n, o, t, r, s, l, a ;计算过程
  ;print, aa ;b, u, p, c, h, n, o, t, r, s, l, a
  ;
  ;---------------------------------------------
  xrdfile = aa
  bgfile = b
  u = u ;refer to the DAT file: 'Um - atomic weight.dat'; rho is g/cm^3
  rho0 = p ;the value of
  ;Scattering factors: 1-C,2-H,3-N,4-0;
  c1 = c & c2 = h & c3 = n & c4 = o; atomic concentration;
  t = t;mm
  R = r;mm
  S = s;rad
  lambda = l; wavelength - Mo=0.701; Cu=1.54
  alpha = a
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
  xrdRead=READ_ASCII(xrdfile)
  ;read backgroudfile.
  bgRead=READ_ASCII(bgfile)
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
;---------------------------------------------
  
  
end

PRO pan
  ;- Create top level base
  tlb = widget_base(column=1,title='PAN-RDF数据处理', $
    tlb_frame_attr=1,XOFFSET=500,YOFFSET=500)
    
  label = widget_label(tlb, value='欢  迎  使  用  ！',/align_center)
  
  ;- Create base to hold everything except buttons
  main = widget_base(tlb, column=1, frame=1)
  ;- Create file widgets
  fbase = widget_base(main, row=1, /base_align_center)
  label = widget_label(fbase, value='XRD数据文件')
  xrd = widget_text(fbase, /editable, xsize=25)
  butt = widget_button(fbase, value='打开文件',UVALUE='xrd')
  
  fbase = widget_base(main, row=1, /base_align_center)
  label = widget_label(fbase, value='XRD背景文件')
  xrdd = widget_text(fbase, /editable, xsize=25)
  butt = widget_button(fbase, value='打开文件',UVALUE='xrdd')
  
  
  
  ;- Create array size widgets
  abase = widget_base(main, row=11, $
    /grid_layout, /base_align_center)
  label = widget_label(abase, value='线吸收系数:')
  u = widget_text(abase, /editable, xsize=8,UVALUE='u')
  
  label = widget_label(abase, value='平均原子数密度（g/cm^3）:')
  p = widget_text(abase, /editable, xsize=8,UVALUE='p')
  
  label = widget_label(abase, value='C原子个数百分比:')
  c = widget_text(abase, /editable, xsize=8,UVALUE='c')
  
  label = widget_label(abase, value='H原子个数百分比:')
  h = widget_text(abase, /editable, xsize=8,UVALUE='h')
  
  label = widget_label(abase, value='N原子个数百分比:')
  n = widget_text(abase, /editable, xsize=8,UVALUE='n')
  
  label = widget_label(abase, value='O原子个数百分比:')
  o= widget_text(abase, /editable, xsize=8,UVALUE='o')
  
  label = widget_label(abase, value='t:')
  t = widget_text(abase, /editable, xsize=8,UVALUE='t')
  
  label = widget_label(abase, value='R:')
  r = widget_text(abase, /editable, xsize=8,UVALUE='r')
  
  label = widget_label(abase, value='S:')
  s = widget_text(abase, /editable, xsize=8,UVALUE='s')
  
  label = widget_label(abase, value='lambda:')
  l = widget_text(abase, /editable, xsize=8,UVALUE='l')
  
  label = widget_label(abase, value='alpha:')
  a = widget_text(abase, /editable, xsize=8,UVALUE='a')
  
  ;- Create ok and cancel buttons
  
  
  bbase = widget_base(tlb, row=1, /align_center)
  buttok = widget_button(bbase, value='计算喽!', xsize=75,UVALUE='run')
  buttres = widget_button(bbase, value='清空', xsize=75,UVALUE='reset')
  ;显示窗口
  widget_control, tlb, /realize
  ;构造参数中介--结构体
  gz={$
    u: u,      $
    p: p,      $
    c: c,      $
    h: h,      $
    n: n,      $
    o: o,      $
    t: t,      $
    r: r,      $
    s: s,      $
    l: l,      $
    a: a,      $
    xrd:xrd,   $
    xrdd:xrdd  $
    
    }
    
  WIDGET_CONTROL,tlb,set_UValue = gz ;将tlb参数信息传递给变量，等待填入数据
  
  ;响应事件
   XMANAGER, 'tlb', tlb, /NO_BLOCK
  
END


;构造XMANAGER函数----------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------


