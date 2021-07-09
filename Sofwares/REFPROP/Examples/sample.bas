Attribute VB_Name = "Module1"
Option Explicit

Public Declare Sub SETUPdll Lib "REFPROP.DLL" (i As Long, ByVal hfld As String, ByVal hfmix As String, ByVal hrf As String, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long, ln3 As Long, ln4 As Long)
Public Declare Sub SETREFdll Lib "REFPROP.DLL" (ByVal hrf As String, ixflag As Long, x0 As Double, h0 As Double, s0 As Double, t0 As Double, p0 As Double, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long)
Public Declare Sub SETMIXdll Lib "REFPROP.DLL" (ByVal hmxnme As String, ByVal hfmix As String, ByVal hrf As String, ncc As Long, ByVal hfile As String, x As Double, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long, ln3 As Long, ln4 As Long, ln5 As Long)
Public Declare Sub SETMODdll Lib "REFPROP.DLL" (i As Long, ByVal htype As String, ByVal hmix As String, ByVal hcomp As String, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long, ln3 As Long, ln4 As Long)
Public Declare Sub SETPATHdll Lib "REFPROP.DLL" (ByVal htype As String, ln As Long)
Public Declare Sub GERG04dll Lib "REFPROP.DLL" (nc As Long, iflag As Long, ierr As Long, ByVal herr As String, ln1 As Long)
Public Declare Sub TPRHOdll Lib "REFPROP.DLL" (t As Double, p As Double, x As Double, j As Long, i As Long, d As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub THERMdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, p As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, hjt As Double)
Public Declare Sub THERM2dll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, p As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, z As Double, hjt As Double, aH As Double, g As Double, kappa As Double, beta As Double, dPdD As Double, d2PdD2 As Double, dPdT As Double, dDdT As Double, dDdP As Double, d2PT2 As Double, d2PdTD As Double, spare3 As Double, spare4 As Double)
Public Declare Sub THERM3dll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, kappa As Double, beta As Double, isenk As Double, kt As Double, betas As Double, bs As Double, kkt As Double, thrott As Double, pi As Double, spht As Double)
Public Declare Sub THERM0dll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, p As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, a As Double, g As Double)
Public Declare Sub ENTROdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, s As Double)
Public Declare Sub ENTHALdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, h As Double)
Public Declare Sub CVCPdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, Cv As Double, Cp As Double)
Public Declare Sub PRESSdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, p As Double)
Public Declare Sub AGdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, a As Double, g As Double)

Public Declare Sub DPDDdll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, dPdD As Double)
Public Declare Sub DPDD2dll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, d2PdD2 As Double)
Public Declare Sub DPDTdll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, dPdT As Double)
Public Declare Sub DDDPdll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, dDdP As Double)
Public Declare Sub DDDTdll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, dDdT As Double)
Public Declare Sub DHD1dll Lib "REFPROP.DLL" (t As Double, rho As Double, x As Double, dHdT_D As Double, dHdT_P As Double, dHdD_T As Double, dHdD_P As Double, dHdP_T As Double, dHdP_D As Double)

Public Declare Sub SATTdll Lib "REFPROP.DLL" (t As Double, x As Double, i As Long, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SATPdll Lib "REFPROP.DLL" (p As Double, x As Double, i As Long, t As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SATDdll Lib "REFPROP.DLL" (d As Double, x As Double, kph As Long, kr As Long, t As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SATHdll Lib "REFPROP.DLL" (h As Double, x As Double, kph As Long, nroot As Long, k1 As Long, t1 As Double, p1 As Double, D1 As Double, k2 As Long, t2 As Double, p2 As Double, d2 As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SATEdll Lib "REFPROP.DLL" (e As Double, x As Double, kph As Long, nroot As Long, k1 As Long, t1 As Double, p1 As Double, D1 As Double, k2 As Long, t2 As Double, p2 As Double, d2 As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SATSdll Lib "REFPROP.DLL" (s As Double, x As Double, kph As Long, nroot As Long, k1 As Long, t1 As Double, p1 As Double, D1 As Double, k2 As Long, t2 As Double, p2 As Double, d2 As Double, k3 As Long, t3 As Double, p3 As Double, d3 As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub CV2PKdll Lib "REFPROP.DLL" (icomp As Long, t As Double, rho As Double, cv2p As Double, csat As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub CSATKdll Lib "REFPROP.DLL" (icomp As Long, t As Double, kph As Long, p As Double, rho As Double, csat As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub DPTSATKdll Lib "REFPROP.DLL" (icomp As Long, t As Double, kph As Long, p As Double, rho As Double, csat As Double, dpdtsat As Double, ierr As Long, ByVal herr As String, ln As Long)

Public Declare Sub TPFLSHdll Lib "REFPROP.DLL" (t As Double, p As Double, x As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub TDFLSHdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub PDFLSHdll Lib "REFPROP.DLL" (p As Double, d As Double, x As Double, t As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub PHFLSHdll Lib "REFPROP.DLL" (p As Double, h As Double, x As Double, t As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub PSFLSHdll Lib "REFPROP.DLL" (p As Double, s As Double, x As Double, t As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub PEFLSHdll Lib "REFPROP.DLL" (p As Double, e As Double, x As Double, t As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub THFLSHdll Lib "REFPROP.DLL" (t As Double, h As Double, x As Double, i As Long, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub TSFLSHdll Lib "REFPROP.DLL" (t As Double, s As Double, x As Double, i As Long, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub TEFLSHdll Lib "REFPROP.DLL" (t As Double, e As Double, x As Double, i As Long, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub DHFLSHdll Lib "REFPROP.DLL" (d As Double, h As Double, x As Double, t As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub DSFLSHdll Lib "REFPROP.DLL" (d As Double, s As Double, x As Double, t As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, h As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub DEFLSHdll Lib "REFPROP.DLL" (d As Double, e As Double, x As Double, t As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub HSFLSHdll Lib "REFPROP.DLL" (h As Double, s As Double, z As Double, t As Double, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, e As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub ESFLSHdll Lib "REFPROP.DLL" (e As Double, s As Double, z As Double, t As Double, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, h As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub CCRITdll Lib "REFPROP.DLL" (t As Double, p As Double, v As Double, x As Double, cs As Double, ts As Double, Ds As Double, ps As Double, ws As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub FPVdll Lib "REFPROP.DLL" (t As Double, d As Double, p As Double, x As Double, f As Double)

Public Declare Sub TQFLSHdll Lib "REFPROP.DLL" (t As Double, q As Double, x As Double, kq As Long, p As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub PQFLSHdll Lib "REFPROP.DLL" (p As Double, q As Double, x As Double, kq As Long, t As Double, d As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub ABFL1dll Lib "REFPROP.DLL" (a As Double, b As Double, x As Double, i As Long, ByVal ab As String, dmin As Double, dmax As Double, t As Double, p As Double, d As Double, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long)
Public Declare Sub ABFL2dll Lib "REFPROP.DLL" (a As Double, b As Double, x As Double, kq As Long, ksat As Long, ByVal ab As String, tbub As Double, tdew As Double, pbub As Double, pdew As Double, Dlbub As Double, Dvdew As Double, ybub As Double, xdew As Double, t As Double, p As Double, Dl As Double, Dv As Double, x As Double, y As Double, q As Double, ierr As Long, ByVal herr As String, ln As Long, ln2 As Long)
Public Declare Sub DBFL2dll Lib "REFPROP.DLL" (d As Double, b As Double, x As Double, i As Long, ByVal ab As String, t As Double, p As Double, Dl As Double, Dv As Double, xliq As Double, xvap As Double, q As Double, ierr As Long, ByVal herr As String, ln As Long, ln2 As Long)

Public Declare Sub CRITPdll Lib "REFPROP.DLL" (x As Double, tc As Double, pc As Double, dc As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub VIRBdll Lib "REFPROP.DLL" (t As Double, x As Double, b As Double)
Public Declare Sub DBDTdll Lib "REFPROP.DLL" (t As Double, x As Double, dbt As Double)
Public Declare Sub VIRCdll Lib "REFPROP.DLL" (t As Double, x As Double, c As Double)
Public Declare Sub TRNPRPdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, eta As Double, tcx As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub FGCTYdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, f As Double)
Public Declare Sub DIELECdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, de As Double)
Public Declare Sub SURFTdll Lib "REFPROP.DLL" (t As Double, d As Double, x As Double, sigma As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SURTENdll Lib "REFPROP.DLL" (t As Double, rhol As Double, rhov As Double, xl As Double, xv As Double, sigma As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub MELTTdll Lib "REFPROP.DLL" (t As Double, x As Double, p As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub MLTH2Odll Lib "REFPROP.DLL" (t As Double, p1 As Double, p2 As Double)
Public Declare Sub MELTPdll Lib "REFPROP.DLL" (p As Double, x As Double, t As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SUBLTdll Lib "REFPROP.DLL" (t As Double, x As Double, p As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub SUBLPdll Lib "REFPROP.DLL" (p As Double, x As Double, t As Double, ierr As Long, ByVal herr As String, ln As Long)

Public Declare Sub WMOLdll Lib "REFPROP.DLL" (x As Double, wm As Double)
Public Declare Sub XMASSdll Lib "REFPROP.DLL" (xmol As Double, xkg As Double, wmix As Double)
Public Declare Sub XMOLEdll Lib "REFPROP.DLL" (xkg As Double, xmol As Double, wmix As Double)
Public Declare Sub QMASSdll Lib "REFPROP.DLL" (qmol As Double, xl As Double, xv As Double, qkg As Double, xlkg As Double, xvkg As Double, wliq As Double, wvap As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub QMOLEdll Lib "REFPROP.DLL" (qkg As Double, xlkg As Double, xvkg As Double, qmol As Double, xl As Double, xv As Double, wliq As Double, wvap As Double, ierr As Long, ByVal herr As String, ln As Long)
Public Declare Sub INFOdll Lib "REFPROP.DLL" (icomp As Long, wmm As Double, ttrp As Double, tnbpt As Double, tc As Double, pc As Double, dc As Double, Zc As Double, acf As Double, dip As Double, Rgas As Double)
Public Declare Sub LIMITXdll Lib "REFPROP.DLL" (ByVal htyp As String, t As Double, d As Double, p As Double, x As Double, tmin As Double, tmax As Double, dmax As Double, pmax As Double, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long)
Public Declare Sub LIMITKdll Lib "REFPROP.DLL" (ByVal htyp As String, icomp As Long, t As Double, d As Double, p As Double, tmin As Double, tmax As Double, dmax As Double, pmax As Double, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long)

Public Declare Sub SETKTVdll Lib "REFPROP.DLL" (icomp As Long, jcomp As Long, ByVal hmodij As String, fij As Double, ByVal hfmix As String, ierr As Long, ByVal herr As String, ln1 As Long, ln2 As Long, ln3 As Long)
Public Declare Sub GETKTVdll Lib "REFPROP.DLL" (icomp As Long, jcomp As Long, ByVal hmodij As String, fij As Double, ByVal hfmix As String, ByVal hfij As String, ByVal hbinp As String, ByVal hmxrul As String, ln1 As Long, ln2 As Long, ln3 As Long, ln4 As Long, ln5 As Long)
Public Declare Sub GETFIJdll Lib "REFPROP.DLL" (ByVal hmodij As String, fij As Double, ByVal hfij As String, ByVal hmxrul As String, ln1 As Long, ln2 As Long, ln3 As Long)
Public Declare Sub PREOSdll Lib "REFPROP.DLL" (i As Long)
Public Declare Sub SETAGAdll Lib "REFPROP.DLL" (ierr As Long, ByVal herr As String, ln1 As Long)


Public Const MaxComps As Integer = 20

Public herr As String * 255, hfmix As String * 255, htyp As String * 3
Public hrf As String * 3
Public htype As String * 3, hmix As String * 3, hcomp As String * 255
Public hfld As String * 10000, hmxnme As String * 255, nc As Long
Public x(1 To MaxComps) As Double, xmm(1 To MaxComps) As Double
Public xliq(1 To MaxComps) As Double, xvap(1 To MaxComps) As Double
Public xv(1 To MaxComps) As Double, rhov As Double, hv As Double
Public ierr As Long, kq As Long, kr As Long, nroot As Long, ix As Long
Public t As Double, p As Double, d As Double, Dl As Double, Dv As Double, q As Double, wm As Double
Public e As Double, h As Double, s As Double, Cv As Double, Cp As Double, w As Double
Public z As Double, hjt As Double, a As Double, g As Double, xkappa As Double, beta As Double, dPdD As Double, d2PdD2 As Double, dPdT As Double, dDdT As Double, dDdP As Double, d2PT2 As Double, d2PdTD As Double, spare3 As Double, spare4 As Double
Public b As Double, c As Double, de As Double, eta As Double, tcx As Double, hvap As Double
Public tmin As Double, tmax As Double, dmax As Double, pmax As Double
Public tc As Double, pc As Double, dc As Double
Public Tref As Double, pref As Double, hRef As Double, sRef As Double
Public f(1 To MaxComps) As Double, sigma As Double
Public k1 As Long, k2 As Long, k3 As Long
Public t1 As Double, p1 As Double, D1 As Double
Public t2 As Double, p2 As Double, d2 As Double
Public t3 As Double, p3 As Double, d3 As Double
Public wmm As Double, ttrp As Double, tnbpt As Double, Zc As Double, acf As Double, dip As Double, Rgas As Double


Sub Example()
  Dim i As Integer
  
  'Pure fluid:
  nc = 1
  hfld = "nitrogen.fld"

  'Mixture:
  'nc = 3
  'hfld = "nitrogen.fld|argon.fld|oxygen.fld|"
  'x(1) = 0.7812
  'x(2) = 0.0092
  'x(3) = 0.2096

  'See fortran for inputs:
  'Call SETMODdll(nc, htype, hmix, hcomp, ierr, herr, 3&, 3&, 255&, 255&)

  hfmix = "hmx.bnc"
  hrf = "DEF"
  Call SETUPdll(nc, hfld, hfmix, hrf, ierr, herr, 10000&, 255&, 3&, 255&)

  'See fortran for inputs:
  'hrf = "OTH"
  'hRef = 1000
  'sRef = 1
  'Tref = 300
  'pref = 1
  'ix = 2
  'Call SETREFdll(hrf, ix, x(1), hRef, sRef, Tref, pref, ierr, herr, 3&, 255&)
  'This is just a check to make sure that h comes back as 1000 J/mol, and s comes back as 1 J/mol-K
  'Call TPFLSHdll(Tref, pref, x(1), d, Dl, Dv, xliq(1), xvap(1), q, e, h, s, cv, cp, w, ierr, herr, 255&)


  t = 100    'K
  p = 1000   'kPa

  'Get fluid info for component 1
  Call INFOdll(1, wmm, ttrp, tnbpt, tc, pc, dc, Zc, acf, dip, Rgas)

  'Get mixture critical parameters:
  Call CRITPdll(x(1), tc, pc, dc, ierr, herr, 255&)

  'Get molecular weight:
  Call WMOLdll(x(1), wm)

  'Get limits:
  htyp = "EOS"
  Call LIMITXdll(htyp, t, d, p, x(1), tmin, tmax, dmax, pmax, ierr, herr, 3&, 255&)

  'Get d from t,p if phase is known:
  'Fourth input: 1 - liquid, 2 - vapor
  'Fourth and fifth inputs must be LONG values!  (include the &)
  Call TPRHOdll(t, p, x(1), 1&, 0&, d, ierr, herr, 255&)

  'Get d from t,p if phase is not known:
  Call TPFLSHdll(t, p, x(1), d, Dl, Dv, xliq(1), xvap(1), q, e, h, s, Cv, Cp, w, ierr, herr, 255&)

  'Calculate properties given T and d:
  Call PRESSdll(t, d, x(1), p)
  Call THERMdll(t, d, x(1), p, e, h, s, Cv, Cp, w, hjt)
  Call THERM2dll(t, d, x(1), p, e, h, s, Cv, Cp, w, z, hjt, a, g, xkappa, beta, dPdD, d2PdD2, dPdT, dDdT, dDdP, d2PT2, d2PdTD, spare3, spare4)
  Call VIRBdll(t, x(1), b)
  Call VIRCdll(t, x(1), c)
  Call DIELECdll(t, d, x(1), de)
  Call FGCTYdll(t, d, x(1), f(1))
  Call SURTENdll(t, d, rhov, x(1), xv(1), sigma, ierr, herr, 255&)
  Call TRNPRPdll(t, d, x(1), eta, tcx, ierr, herr, 255&)

  'Saturation routines:
  'Third input for Satt and Satp is:  1 - bubble point, 2 - dew point
  Call SATTdll(t, x(1), 1&, p, Dl, Dv, xliq(1), xvap(1), ierr, herr, 255&)
  Call SATPdll(p, x(1), 1&, t, Dl, Dv, xliq(1), xvap(1), ierr, herr, 255&)
  'Call SATDdll(d, x(1), 1&, kr, t, p, Dl, Dv, xliq(1), xvap(1), ierr, herr, 255&)
  'Call SATHdll(h, x(1), 0&, nroot, k1, t, p1, d1, k2, t2, p2, d2, ierr, herr, 255&)
  'Call SATSdll(s, x(1), 0&, nroot, k1, t1, p1, d1, k2, t2, p2, d2, k3, t3, p3, d3, ierr, herr, 255&)

  'Routines to calculate properties when T and d are not both known:
  'Call DEFLSHdll(d, e, x(1), t, p, Dl, Dv, xliq(1), xvap(1), q, h, s, cv, cp, w, ierr, herr, 255&)
  'Call DHFLSHdll(d, h, x(1), t, p, Dl, Dv, xliq(1), xvap(1), q, e, s, cv, cp, w, ierr, herr, 255&)
  'Call DSFLSHdll(d, s, x(1), t, p, Dl, Dv, xliq(1), xvap(1), q, e, h, cv, cp, w, ierr, herr, 255&)
  'Call PDFLSHdll(p, d, x(1), t, Dl, Dv, xliq(1), xvap(1), q, e, h, s, cv, cp, w, ierr, herr, 255&)
  'Call PEFLSHdll(p, e, x(1), t, d, Dl, Dv, xliq(1), xvap(1), q, e, h, cv, cp, w, ierr, herr, 255&)
  'Call PHFLSHdll(p, h, x(1), t, d, Dl, Dv, xliq(1), xvap(1), q, e, s, cv, cp, w, ierr, herr, 255&)
  'Call PQFLSHdll(p, q, x(1), 1&, t, d, Dl, Dv, xliq(1), xvap(1), e, h, s, cv, cp, w, ierr, herr, 255&)
  'Call PSFLSHdll(p, s, x(1), t, d, Dl, Dv, xliq(1), xvap(1), q, e, h, cv, cp, w, ierr, herr, 255&)
  'Call TDFLSHdll(t, d, x(1), p, Dl, Dv, xliq(1), xvap(1), q, e, h, s, cv, cp, w, ierr, herr, 255&)
  'Call TEFLSHdll(t, e, x(1), 1&, p, d, Dl, Dv, xliq(1), xvap(1), q, h, s, cv, cp, w, ierr, herr, 255&)
  'Call THFLSHdll(t, h, x(1), 1&, p, d, Dl, Dv, xliq(1), xvap(1), q, e, s, cv, cp, w, ierr, herr, 255&)
  'Call TQFLSHdll(t, q, x(1), 1&, p, d, Dl, Dv, xliq(1), xvap(1), e, h, s, cv, cp, w, ierr, herr, 255&)
  'Call TSFLSHdll(t, s, x(1), 1&, p, d, Dl, Dv, xliq(1), xvap(1), q, e, h, cv, cp, w, ierr, herr, 255&)

  'Melting and sublimation routines:
  'Call MELTPdll(p, x(1), t, ierr, herr, 255&)
  'Call MELTTdll(t, x(1), p, ierr, herr, 255&)
  'Call SUBLTdll(t, x(1), p, ierr, herr, 255&)
  'Call SUBLPdll(p, x(1), t, ierr, herr, 255&)

  'Convert mole fraction to/from mass fraction:
  'Call XMASSdll(x(1), xmm(1), wm)
  'Call XMOLEdll(x(1), xmm(1), wm)


  'Example setup for R404A:
  nc = 3
  hfld = "C:\Program Files\REFPROP\fluids\R125.FLD|C:\Program Files\REFPROP\fluids\R134A.FLD|C:\Program Files\REFPROP\fluids\R143A.FLD|"
  x(1) = 0.357816784026318
  x(2) = 3.82639950410712E-02
  x(3) = 0.603919220932611
  hfmix = "C:\Program Files\REFPROP\fluids\hmx.bnc"
  hrf = "DEF"
  Call SETUPdll(nc, hfld, hfmix, hrf, ierr, herr, 10000&, 255&, 3&, 255&)
  t = 300     'K
  p = 10000   'kPa
  Call TPFLSHdll(t, p, x(1), d, Dl, Dv, xliq(1), xvap(1), q, e, h, s, Cv, Cp, w, ierr, herr, 255&)
  'Answer for d should be 11.371439 mol/dm^3 and for w should be 488.5738 m/s

  'In version 8.0, R404A can also be setup with a call to SETMIX.
  'In this situation, the values in the x array are returned by SETMIX
  'and do not need to be set initially.  Only hmxnme, hfmix, and hrf are inputs.
  hmxnme = "C:\Program Files\REFPROP\mixtures\R404A.mix"
  hfmix = "C:\Program Files\REFPROP\fluids\hmx.bnc"
  hrf = "DEF"
  Call SETMIXdll(hmxnme, hfmix, hrf, nc, hfld, x(1), ierr, herr, 255&, 255&, 3&, 10000&, 255&)
  t = 300     'K
  p = 10000   'kPa
  Call TPFLSHdll(t, p, x(1), d, Dl, Dv, xliq(1), xvap(1), q, e, h, s, Cv, Cp, w, ierr, herr, 255&)
End Sub
