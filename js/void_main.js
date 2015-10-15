// Void Evolution in Pure Aluminum
// Javascript functions 
// Tianyu Liu

var $ = function(id) {
  return document.getElementById(id);
}

function Evolute() {
  outEv = [];
  T = $("tT").value, Ds = calcDs()
  var tEv = $('tEv').value, tEv0 = tEv
  var gph_opt = $("gphEv_opt").value
  var outj1 = [], dn = 0, fMax = 0, f = 0, dt = 0, tMin = 0, i0 = 1, Ii = [], Iu =
      [], Iui = [ 1, -1 ]

  for (var k = 1; k <= 14; k++)
    Iu[k] = 0;

  Ii = getAllI()
  outEv[0] = cpEv(EvSubStep(Ii))

  logEv +=
      "Step #: Vacancies, Total Surface Energy; 14 Face distances (atomic layer),  D{002}/D{111}\n\n"
  logEv += "   Sub-step. This face to act on: ΔE/ΔN; D{002}/D{111}\n"
  logEv += " * (This sub-step has the largest ΔE/ΔN so far)\n"
  logEv +=
      " *^(Max ΔE/ΔN too, and this face hasn't been updated for too long)\n\n"
  logEv +=
      " $ Selected sub-step, ΔE/ΔN, Steps since each face's opposite been grown (cut)\n\n"

  var j = -1
  do {
    j++
    if (tEv0 == 0) tEv++

    logEv +=
        "\nStep " + j + ": " + outEv[j].n + ", " + outEv[j].E.toExponential(4)
            + ";  "
    for (var k = 1; k <= 14; k++) {
      logEv += outEv[j].I[k] + ","
    }
    logEv += " " + outEv[j].I[0].toFixed(5) + "\n"

    if (gph_opt == 0) outEv[j].P = genCoor(outEv[j].I);
    fMax = 0, i0 = 1

    for (var i = 1; i <= 14; i++) {
      Ii = outEv[j].I.slice(), Ii[i]--, Ii = calcI(Ii)
      outj1[i] = EvSubStep(Ii)
      dn = outEv[j].n - outj1[i].n
      f = (outEv[j].E - outj1[i].E) / dn

      if (f > fMax) {
        i0 = i
        fMax = f
        logEv += "\n   * "
      } else if (f == fMax && Iu[i] > Iu[i0]) {
        // avoid acting only on some faces
        i0 = i
        fMax = f
        logEv += "\n   *^"
      } else logEv += "\n     ";

      logEv +=
          j + "-" + fmtStr(i + ".", 3) + FI[i] + ": " + f.toExponential(3)
              + "; "
      logEv += Ii[0].toFixed(5)
    }

    outEv[j].dn = outEv[j].n - outj1[i0].n, outEv[j].f = fMax, outEv[j].i = i0
    outEv[j].inR = inRatio(outEv[j])
    outEv[j].dt = calcTime(outEv[j])
    outEv[j + 1] = cpEv(outj1[i0])

    // Iu records how long since the opposite of each face last update
    for (var k = 1; k <= 14; k++)
      if (outEv[j].I[k] != outEv[j + 1].I[k]) {
        Iu[k - Iui[k % 2]]++;
        Iu[k]--;
      }
    ;

    logEv +=
        "\n\n   $ " + j + "-" + fmtStr(i0 + ".", 3) + FI[i0] + ": "
            + fMax.toExponential(3) + ";  "
    for (var k = 1; k <= 14; k++)
      logEv += Iu[k] + ",";
  } while (outEv[j + 1].n > 0 && outEv[j + 1].E > 0 && outEv[j].dt > 0
      && j < tEv)

  outEv.pop()
  tEv = j

  textEvSmy()
  outputField(logEv, "field_EvLog")
  outputField(smyEv, "field_EvSumm")
  $("iEv_opt").options.length = 0
  for (var i = 0; i <= tEv; i++)
    $("iEv_opt").options[i] = new Option(i, i);
  iEv = tEv;
  dispEv();
}

function EvSubStep(I) {
  var n = calcAtom(I)
  var S = calcS(I)

  return {
    n : n,
    E : S[0],
    I : I,
    S : S
  }
}

function cpEv(out) {
  return {
    n : out.n,
    E : out.E,
    I : out.I,
    S : out.S,
    i : 0,
    dn : 0,
    f : 0,
    inR : [],
    dt : 0,
    P : []
  }
}

function inRatio(out) {
  var i = out.i, I1 = out.I.slice(), inR = []

  for (var k = 9; k <= 14; k++)
    I1[k] = Math.round(I1[k] * t21);

  for (var k = 1; k <= 14; k++) {
    if (Math.abs(I1[k] - I1[i]) < 5) inR[k] = "Y";
    else inR[k] = "N";
  }

  return inR
}

function calcTime(out) {
  // Calculate time
  var S_D = 0, I = out.I, S = out.S, dn = out.dn, f = out.f, inR = out.inR

  for (var k = 1; k <= 8; k++) {
    if (inR[k] == "Y") S_D += S[k] / I[k] * sqrt3;
  }
  for (var k = 9; k <= 14; k++) {
    if (inR[k] == "Y") S_D += S[k] / I[k] * 2;
  }

  var dt = dn / (S_D * a / 10 * Ds / omega / xi * (Math.exp(f / kB / T) - 1))
  return dt
}

function textEvSmy() {
  var i = 0, dt = 0, time = 0, inR = []

  smyEv = "", time = 0, Ds = calcDs(), T = $("tT").value, Q = $("tQ").value,
      D0 = $("tD0").value
  smyEv += "Ds\tT\tQ\tD0\n"
  smyEv += "nm²/s\tK\teV\tnm²/s\n"
  smyEv += Ds.toExponential(3) + "\t" + T + "\t"
  if (setDs() == 0) smyEv += Q + "\t" + D0
  else smyEv += "-\t-"

  smyEv += "\n\nStep\tSelected\t"
  smyEv += "Vacancies\tSurface E\tΔN\tΔE/ΔN\tΔt\tTime\tD{111}\tD{002}\tRatio"
  for (k = 1; k <= 14; k++)
    smyEv += "\t" + k + ". I" + FI[k];
  for (k = 1; k <= 14; k++)
    smyEv += "\t" + k + ". S" + FI[k];
  smyEv += "\n \t \t \teV\t\teV/vacancy\tmin\tmin\tnm\tnm\t "
  for (k = 1; k <= 14; k++)
    smyEv += "\tAtom I.";
  for (k = 1; k <= 14; k++)
    smyEv += "\tnm²";
  for (k = 1; k <= 14; k++)
    smyEv += "\t" + k + ". " + FI[k];

  for (var j = 0; j <= outEv.length - 1; j++) {
    i = outEv[j].i, dt = outEv[j].dt, time += dt, I = outEv[j].I, inR =
        outEv[j].inR

    smyEv += "\n" + j + "\t" + i + ". '" + FI[i] + "'\t"
    smyEv +=
        outEv[j].n + "\t" + outEv[j].E.toExponential(15) + "\t" + outEv[j].dn
            + "\t"
    smyEv +=
        outEv[j].f.toExponential(15) + "\t" + (dt / 60).toExponential(15)
            + "\t"
    smyEv +=
        (time / 60).toExponential(15) + "\t"
            + (outEv[j].I[15] / sqrt3 * a / 5).toExponential(15) + "\t"
    smyEv +=
        (outEv[j].I[16] * a / 10).toExponential(15) + "\t"
            + outEv[j].I[0].toExponential(15)
    for (var k = 1; k <= 14; k++)
      smyEv += "\t" + outEv[j].I[k];
    for (var k = 1; k <= 14; k++)
      smyEv += "\t" + (outEv[j].S[k] * as).toFixed(5);
    for (var k = 1; k <= 14; k++)
      smyEv += "\t" + inR[k];
  }
}

function calcMaxDif(I) {
  var I1 = I.slice(), max, min

  for (var k = 9; k <= 14; k++)
    I1[k] = Math.round(I[k] * t21);

  max = min = I1[1]
  for (var k = 2; k <= 14; k++) {
    if (I1[k] > max) max = I1[k];
    else if (I1[k] < min) min = I1[k];
  }

  dmax = max - min

  return dmax
}

function dispEv() {
  var tEv = outEv.length - 1

  if (0 <= iEv && iEv <= tEv) {
    var I = outEv[iEv].I;
    P = outEv[iEv].P;

    outputField("", "field_AtomCd")

    dispSpec(outEv[iEv].n, I, outEv[iEv].S)
    if ($("gphEv_opt").value == 1) P = genCoor(I);

    if ($('atomCd_chk').checked == true) dispCoor(P);
    draw3D(P);

    $('iEv_opt').options[iEv].selected = 1
  } else if (iEv < 0) iEv = 0;
  else iEv = tEv
}

function Display() {
  var rt = [ performance.now() ];

  var I = getAllI();
  rt[1] = performance.now();
  dispSpec(calcAtom(I), I, calcS(I));
  rt[2] = performance.now();
  P = genCoor(I);
  rt[3] = performance.now();
  if ($('atomCd_chk').checked == true) dispCoor(P);
  rt[4] = performance.now();
  draw3D(P);
  rt[5] = performance.now();
  rt[6] = (rt[5] - rt[0]).toFixed(3);

  for (var i = 0; i < 5; i++)
    rt[i] = (rt[i + 1] - rt[i]).toFixed(3);
  console.log("Display finished in ", +rt[6], " ms: getAllI ", +rt[0],
      ", dispSpec ", +rt[1], ", genCoor ", +rt[2], ", dispCoor ", +rt[3],
      ", draw3D ", +rt[4]);
}

function dispCoor(P) {
  var Cd = P[0], CdN = P[1], mL = P[2];
  var i = j = 0, n_Vd = 0, n_Mtx = CdN[0];
  textVd = textMtxVd = textMtx = "";

  for (i = 0; i < n_Mtx; i++)
    if (Cd[0][i][0] !== mL && Cd[0][i][1] !== mL && Cd[0][i][2] !== mL)
      textMtx += (i + 1) + " 1 " + fmtcdNum(Cd[0][i], ".5") + "\n";

  textMtxVd = textMtx;
  var texti = "";

  for (j = 1; j < 16; j++)
    for (i = 0; i < CdN[j]; i++) {
      texti = " 1 " + fmtcdNum(Cd[j][i], ".5") + "\n";
      textVd += ++n_Vd + texti;
      textMtxVd += (n_Vd + n_Mtx) + texti
    }

  mL = (mL * a).toFixed(5);

  var textHead = "1 atom types\n";
  textHead += mL + " " + mL + " xlo xhi\n";
  textHead += mL + " " + mL + " ylo yhi\n";
  textHead += mL + " " + mL + " zlo zhi\n";
  textHead += "\nAtoms\n\n";

  textMtx = textHead + textMtx;
  textMtx = n_Mtx + " atoms\n" + textMtx;
  textMtx = "# Matrix (with Void) coordinates for LAMMPS\n" + textMtx;
  textVd = textHead + textVd;
  textVd = n_Vd + " atoms\n" + textVd;
  textVd = "# Void coordinates for LAMMPS\n" + textVd;
  textMtxVd = textHead + textMtxVd;
  textMtxVd = (n_Mtx + n_Vd) + " atoms\n" + textMtxVd;
  textMtxVd = "# Matrix without Void coordinates for LAMMPS\n" + textMtxVd;

  outputField(textMtx, "field_AtomCd")
}

function dispSpec(n, I, S) {
  var NL = 2, Sa, Ia, E

  for (var i = 1; i <= 16; i++) {
    Sa = S[i] * as

    if (i <= 8 || i == 15) {
      Ia = I[i] * a / sqrt3
      E = Sa * 0.06552
    } else {
      Ia = I[i] * a / 2;
      E = Sa * 0.06864
    }

    $("dpIr_" + i).innerHTML = I[i];
    $("dpIa_" + i).innerHTML = Ia.toFixed(NL);
    $("dpSa_" + i).innerHTML = Sa.toExponential(NL);
    $("dpE_" + i).innerHTML = E.toExponential(3)
  }

  $("dpIr_15").innerHTML = I[15].toFixed(NL);
  $("dpIr_16").innerHTML = I[16].toFixed(NL);

  $("dpIr_17").innerHTML = I[0].toFixed(5);
  $("dpSa_17").innerHTML = (S[17] * as).toExponential(NL);
  $("dpE_17").innerHTML = S[0].toExponential(3);
  $("dpN").innerHTML = n;
}

function genCoor(I) {
  var x, y, z, matrixL = Number($("matrixL").value);
  var Cd = [], Cdi = 0, n = 0, CdN = [];
  var mL = Math.max(I[9], I[10], I[11], I[12], I[13], I[14]);

  for (var i = 0; i < 16; i++)
    Cd[i] = [], CdN[i] = 0;

  if (mL >= matrixL) mL = ++mL / 2;
  else mL = matrixL / 2;

  for (x = -mL; x <= mL; x += .5) {
    for (y = -mL; y <= mL; y += .5) {
      for (z = -mL; z <= mL; z += .5) {
        if ((x % 1 + y % 1 + z % 1) % 1 === 0) {
          Cdi = 0;
          if (I[9] / 2 >= x && I[10] / 2 >= -x && I[11] / 2 >= y
              && I[12] / 2 >= -y && I[13] / 2 >= z && I[14] / 2 >= -z
              && I[1] >= +x + y + z && I[2] >= -x - y - z && I[3] >= -x + y + z
              && I[4] >= +x - y - z && I[5] >= +x - y + z && I[6] >= -x + y - z
              && I[7] >= +x + y - z && I[8] >= -x - y + z) {
            Cdi++;
            if (I[1] === x + y + z) Cdi++;
            if (I[2] === -x - y - z) Cdi++;
            if (I[3] === -x + y + z) Cdi++;
            if (I[4] === x - y - z) Cdi++;
            if (I[5] === x - y + z) Cdi++;
            if (I[6] === -x + y - z) Cdi++;
            if (I[7] === x + y - z) Cdi++;
            if (I[8] === -x - y + z) Cdi++;
            if (I[9] / 2 === x) Cdi++;
            if (I[10] / 2 === -x) Cdi++;
            if (I[11] / 2 === y) Cdi++;
            if (I[12] / 2 === -y) Cdi++;
            if (I[13] / 2 === z) Cdi++;
            if (I[14] / 2 === -z) Cdi++;
            // console.log(Cdi, x, y, z)
          }

          Cd[Cdi][CdN[Cdi]++] = [ x, y, z ];
        }
      }
    }
  }

  return [ Cd, CdN, mL ];
}

function setAllI(o) {
  if (o == 0) for (i = 1; i <= 8; i++)
    $("tI_" + i).value = $("tI1_all").value
  else for (i = 9; i <= 14; i++)
    $("tI_" + i).value = $("tI2_all").value
}

function getAllI() {
  var I = []
  for (var i = 1; i <= 14; i++)
    I[i] = parseInt($("tI_" + i).value)

  return calcI(I)
}

function calcI(I) {
  var I1 = [], la1 = la2 = 0

  I1[1] = Math.min(I002toI111(9, 11, 13, I), I[1])
  I1[2] = Math.min(I002toI111(10, 12, 14, I), I[2])
  I1[3] = Math.min(I002toI111(10, 11, 13, I), I[3])
  I1[4] = Math.min(I002toI111(9, 12, 14, I), I[4])
  I1[5] = Math.min(I002toI111(9, 12, 13, I), I[5])
  I1[6] = Math.min(I002toI111(10, 11, 14, I), I[6])
  I1[7] = Math.min(I002toI111(9, 11, 14, I), I[7])
  I1[8] = Math.min(I002toI111(10, 12, 13, I), I[8])
  I1[9] = Math.min(I[1] + I[4], I[5] + I[7], I[9])
  I1[10] = Math.min(I[2] + I[3], I[6] + I[8], I[10])
  I1[11] = Math.min(I[1] + I[6], I[3] + I[7], I[11])
  I1[12] = Math.min(I[2] + I[5], I[4] + I[8], I[12])
  I1[13] = Math.min(I[1] + I[8], I[3] + I[5], I[13])
  I1[14] = Math.min(I[2] + I[7], I[4] + I[6], I[14])

  for (var i = 1; i <= 8; i++)
    la1 += I1[i];
  for (var i = 9; i <= 14; i++)
    la2 += I1[i];

  I1[15] = la1 / 8
  I1[16] = la2 / 6
  I1[0] = (I1[16] / 2) / (I1[15] / sqrt3)

  return I1
}

function calcD(I) {
  var D = []

  for (var i = 1; i <= 8; i++)
    D[i] = I[i] * a_sqrt3;
  for (var i = 9; i <= 14; i++)
    D[i] = I[i] * a_2;

  D[15] = I[15] * a_sqrt3
  D[16] = I[16] * a_2

  return D
}

function I002toI111(a, b, c, I) {
  var x = I[a] + I[b] + I[c];
  return Math.floor(x / 2)
}

function calcS(I) {
  var MS1 = [], MS2 = [], S = [], b = sqrt3 / 8

  MS1[1] = I[1] - I[3] - I[5] - I[7]
  MS1[2] = I[2] - I[4] - I[6] - I[8]
  MS1[3] = I[3] - I[1] - I[6] - I[8]
  MS1[4] = I[4] - I[7] - I[5] - I[2]
  MS1[5] = I[5] - I[1] - I[4] - I[8]
  MS1[6] = I[6] - I[7] - I[3] - I[2]
  MS1[7] = I[7] - I[1] - I[4] - I[6]
  MS1[8] = I[8] - I[3] - I[5] - I[2]

  MS2[1] = I[9] - I[1] - I[4]
  MS2[2] = I[9] - I[7] - I[5]
  MS2[3] = I[10] - I[3] - I[2]
  MS2[4] = I[10] - I[6] - I[8]
  MS2[5] = I[11] - I[1] - I[6]
  MS2[6] = I[11] - I[7] - I[3]
  MS2[7] = I[12] - I[5] - I[2]
  MS2[8] = I[12] - I[4] - I[8]
  MS2[9] = I[13] - I[1] - I[8]
  MS2[10] = I[13] - I[3] - I[5]
  MS2[11] = I[14] - I[2] - I[7]
  MS2[12] = I[14] - I[4] - I[6]

  S[9] = MS2[1] * MS2[2] / 2
  S[10] = MS2[3] * MS2[4] / 2
  S[11] = MS2[5] * MS2[6] / 2
  S[12] = MS2[7] * MS2[8] / 2
  S[13] = MS2[9] * MS2[10] / 2
  S[14] = MS2[11] * MS2[12] / 2

  for (var i = 1; i <= 8; i++)
    MS1[i] *= MS1[i];
  for (var i = 1; i <= 12; i++)
    MS2[i] *= MS2[i];

  S[1] = (MS1[1] - MS2[2] - MS2[6] - MS2[10]) * b
  S[2] = (MS1[2] - MS2[4] - MS2[8] - MS2[12]) * b
  S[3] = (MS1[3] - MS2[4] - MS2[5] - MS2[9]) * b
  S[4] = (MS1[4] - MS2[2] - MS2[7] - MS2[11]) * b
  S[5] = (MS1[5] - MS2[1] - MS2[8] - MS2[9]) * b
  S[6] = (MS1[6] - MS2[3] - MS2[6] - MS2[11]) * b
  S[7] = (MS1[7] - MS2[1] - MS2[5] - MS2[12]) * b
  S[8] = (MS1[8] - MS2[3] - MS2[7] - MS2[10]) * b

  for (var i = 1; i <= 14; i++)
    if (S[i] < 0) S[i] = 0;

  S[15] = 0, S[16] = 0

  for (var i = 1; i <= 8; i++)
    S[15] += S[i];
  for (var i = 9; i <= 14; i++)
    S[16] += S[i];

  S[17] = S[15] + S[16]
  S[0] = (S[15] * Epa1 + S[16] * Epa2) * as

  return S
}

function calcArea(S) {
  var SA = []

  for (var i = 1; i <= 17; i++)
    SA[i] = S[i] * as;

  return SA
}

function calcAtom(I) {
  var n = 0, n0 = 0, np = 1, ns = 0, n1 = 0, n2 = 0, oe = []
  n0 = (I[9] + I[10] + 1) * (I[11] + I[12] + 1) * (I[13] + I[14] + 1) / 2
  for (var i = 9; i <= 14; i++) {
    oe[i] = I[i] % 2;
    if (oe[i] == 1) np = 0
  }
  if (np == 0) {
    for (var i = 9; i <= 14; i++)
      ns += oe[i];
    if (ns == 4) {
      if (oe[9] + oe[10] == 0 || oe[11] + oe[12] == 0 || oe[13] + oe[14] == 0)
        np = 1
    }
  }

  n0 = Math.floor(n0) + np

  n1 += calc111Atom(I[1], I[9], I[11], I[13])
  n1 += calc111Atom(I[2], I[10], I[12], I[14])
  n1 += calc111Atom(I[3], I[10], I[11], I[13])
  n1 += calc111Atom(I[4], I[9], I[12], I[14])
  n1 += calc111Atom(I[5], I[9], I[12], I[13])
  n1 += calc111Atom(I[6], I[10], I[11], I[14])
  n1 += calc111Atom(I[7], I[9], I[11], I[14])
  n1 += calc111Atom(I[8], I[10], I[12], I[13])

  n2 += calc111OverlapAtom(I[1], I[7], I[9], I[11])
  n2 += calc111OverlapAtom(I[1], I[5], I[9], I[13])
  n2 += calc111OverlapAtom(I[1], I[3], I[11], I[13])
  n2 += calc111OverlapAtom(I[7], I[4], I[9], I[14])
  n2 += calc111OverlapAtom(I[7], I[6], I[11], I[14])
  n2 += calc111OverlapAtom(I[5], I[4], I[9], I[12])
  n2 += calc111OverlapAtom(I[5], I[8], I[12], I[13])
  n2 += calc111OverlapAtom(I[4], I[2], I[12], I[14])
  n2 += calc111OverlapAtom(I[3], I[6], I[10], I[11])
  n2 += calc111OverlapAtom(I[3], I[8], I[10], I[13])
  n2 += calc111OverlapAtom(I[6], I[2], I[10], I[14])
  n2 += calc111OverlapAtom(I[8], I[2], I[10], I[12])

  n = n0 - n1 + n2

  return n
}

function calc111OverlapAtom(x, y, a, b) {
  var n0 = calc111Atom(x, a, b, x - y - 1) + calc111Atom(y, a, b, y - x - 1)
  var l = a + b - x - y - 1
  if (l <= 0) return n0;
  return n0 + (l * (2 + l) + l % 2) / 4
}

function calc111Atom(x, a, b, c) {
  var l = Math.floor((a + b + c) / 2) - x
  if (l <= 0) return 0;

  var n = l * (l + 1) / 6

  if ((a + b + c) % 2 == 0) n *= 4 * l - 1
  else n *= 4 * l + 5

  return n
}

function calcDs() {
  if (1 == setDs()) return 1

  var T = $("tT").value, Q = $("tQ").value, D0 = $("tD0").value, Ds =
      D0 * Math.exp(-Q / kB / T)
  $('dDs').innerHTML = Ds.toExponential(3)

  return Ds
}

function setDs() {
  var o = $('setDs_opt').checked

  $('tD0').disabled = o
  $('tQ').disabled = o
  if (o == 1) $('dDs').style = "color:gray"
  else $('dDs').style = "color:auto"

  return o
}

function outputField(text, field) {
  if (text.length < 1e6) $(field).value = text;
  else $(field).value = "Content is too large to display, please download.";
}

function atomButton(filename, text) {
  if ($("dl_chk").checked == true) download(filename, text);
  else $('field_AtomCd').value = text;
}

function download(filename, text) {
  var blob = new Blob([ text ], {
    type : "text/plain;charset=utf-8"
  });
  saveAs(blob, filename + ".txt");
}

function fmtcdNum(Cdi, fmt) {
  var fmtcd = ""
  for (var i = 0; i < 3; i++)
    fmtcd += fmtNum(Cdi[i] * a, fmt) + " ";
  return fmtcd
}

function fmtNum(num, fmt) {
  var fmt = String(fmt), m = fmt.split(".")[0]
  var fmt1 = String(num), m1 = fmt1.split(".")[0]
  num = num.toFixed(fmt.split(".")[1])
  if (m1.length < m) num = Array(m - m1.length + 1).join(" ") + num
  return num
}

function fmtStr(str, fmt) {
  if (str.length < fmt) str += Array(fmt - str.length + 1).join(" ")
  return str
}