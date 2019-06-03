#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

struct scr_resultA{
  int nCyto, nColpo, result;
};

struct scr_resultB{
  int nCyto, nHPV, nColpo, nBiops, nCombos, result;
};

struct scr_resultC{
  int nCyto, nHPV, nColpo, nBiops, nCombos, result;
};

struct scr_resultD{
  int nCyto, nHPV, nColpo, nBiops, nCombos, result;
};

int cytologyTEST(SEXP SENSI, int location, double probscyto[3], int pathPROTO)
{
  int detectedW=0, res=0;
  
  detectedW  = rbinom(1, *REAL(VECTOR_ELT(SENSI, location)));
  if (detectedW == 0)
  {
    detectedW = rbinom(1, 0.039); /* Average ASCUS+LSIL+HSIL */
  }
  if (pathPROTO == 1)
  {
    detectedW = 1;
  }
  if (detectedW == 0)
  {
    res = 0;
  }
  if (detectedW == 1)
  {
    if ((location == 0) | (location == 1))
    {
      res = 1; /* ASCUS */
    }
    if ((location == 2) | (location == 3) | (location == 4))
    {
      res = 2; /* LSIL */
    }
    if ((location == 5) | (location == 6) | (location == 7) | (location == 8))
    {
      res = 3; /* HSIL */
    }
  }
  return res;
}

int colpoTEST(SEXP SENSI, int location)
{
  int res=0;
  
  res = rbinom(1, *REAL(VECTOR_ELT(SENSI, location)));
  return res;
}

struct scr_resultA scrSchema_A(SEXP SCREENSENSI, SEXP SCREENSENSI2, SEXP COLPOSENSI, int location, int pos,
                               double probscyto[3])
{
  int detectedW=0, resCyto2=0, resColpo=0;
  struct scr_resultA res={.nCyto=0, .nColpo=0, .result=0};
  
  if (pos == 1)
  {
    goto CYTO2;
  }
  detectedW = cytologyTEST(SCREENSENSI, location, probscyto, 0);
  res.nCyto = res.nCyto + 1;
  if (detectedW == 1) /* ASCUS */
  {
    CYTO2:resCyto2 = cytologyTEST(SCREENSENSI2, location, probscyto, 0);
    res.nCyto = res.nCyto + 1;
  }
  if ((detectedW == 2) | (detectedW == 3) | (resCyto2 > 0)) /* Neg */
  {
    resColpo = colpoTEST(COLPOSENSI, location);
    res.nColpo = res.nColpo + 1;
    if (resColpo == 0) /* Neg or CIN1 */
    {
      res.result = 2;
    }
    if (resColpo == 1) /* CIN2+ */
    {
      res.result = 3;
    }
  }
  return res;
}

int hpvTEST(SEXP SENSI, int location)
{
  int res=0;
  
  res = rbinom(1, *REAL(VECTOR_ELT(SENSI, location)));
  return res;
}

int colpoPROTO(SEXP SCREENSENSI, SEXP HPVTESTSENSI, SEXP COLPOSENSI, int location, int ColpoBefore)
{
  int detectedW=0, res=0;
  
  detectedW = colpoTEST(COLPOSENSI, location);
  if (detectedW == 0) /* Neg */
  {
    res = 99;
  }
  if (detectedW == 1) /* CIN1+ */
  {
    res = 999;
  }
  return res;
}

int cytohpvPROTO(SEXP SCREENSENSI, SEXP HPVTESTSENSI, SEXP COLPOSENSI, int location, int ColpoBefore,
                 double probscyto[3])
{
  int detectedW=0, detectedW2=0, res=0;
  
  detectedW  = cytologyTEST(SCREENSENSI, location, probscyto, 0);
  detectedW2 = hpvTEST(HPVTESTSENSI, location);
  if ((detectedW == 0) & (detectedW2 == 0)) /* Neg */
  {
    res = 88;
  }
  if ((detectedW != 0) | (detectedW2 == 1)) /* any Pos */
  {
    if (ColpoBefore == 0) /* No colpo before */
    {
      res = 888;
    }
    if (ColpoBefore != 0) /* colpo before */
    {
      res = 8888;
    }
  }
  return res;
}

int biopsPROTO(SEXP SCREENSENSI, SEXP HPVTESTSENSI, SEXP COLPOSENSI, SEXP BIOPSENSI, int location,
               int ColpoBefore)
{
  int detectedW=0, res=0;
  
  detectedW = rbinom(1, *REAL(VECTOR_ELT(BIOPSENSI, location)));
  if (detectedW == 0) /* CIN1 */
  {
    res = 77;
  }
  if (detectedW == 1) /* CIN2-3+ */
  {
    res = 777;
  }
  return res;
}

struct scr_resultB scrSchema_B(SEXP SCREENSENSI, SEXP SCREENSENSI2, SEXP SCREENSENSI3, SEXP COLPOSENSI,
                               SEXP BIOPSENSI, SEXP HPVTESTSENSI, int location, double probscyto[3], int pathPROTO,
                               int cytohpv_TOCA, int cyto_TOCA, int cytoType)
{
  int resCyto1=0, resHpvTest=0, ColpoBefore=0;
  struct scr_resultB res={.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  
  if ((cyto_TOCA == 1) | (cytohpv_TOCA == 1))
  {
    if (cyto_TOCA == 1)
    {
      resCyto1  = cytologyTEST(SCREENSENSI, location, probscyto, 0);
      res.nCyto = res.nCyto + 1;
      if (resCyto1==0)
      {
        res.result = 2;
      }
      if (resCyto1 != 0)
      {
        pathPROTO = 1;
        goto CYTO1;
      }
    }
    if (cytohpv_TOCA == 1)
    {
      res.result = cytohpvPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore, probscyto);
      res.nCyto  = res.nCyto + 1;
      res.nHPV   = res.nHPV + 1;
      if (cytoType == 1)
      {
        res.nCombos = res.nCombos + 1;
      }
      if (res.result == 88)
      {
        res.result = 0;
      }
      if (res.result == 8888)
      {
        res.result = colpoPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
  }else{
    CYTO1:resCyto1 = cytologyTEST(SCREENSENSI, location, probscyto, pathPROTO);
    if (pathPROTO == 0)
    {
      res.nCyto = res.nCyto + 1;
    }
    pathPROTO = 0;
    if (resCyto1 == 1) /* ASCUS */
    {
      resHpvTest = hpvTEST(HPVTESTSENSI, location);
      res.nHPV = res.nHPV + 1;
      if (cytoType == 1)
      {
        res.nCombos = res.nCombos + 1;
      }
      if (resHpvTest == 0) /* HPV- or LR */
      {
        res.result = 1;
      }
      if (resHpvTest == 1) /* HR HPV */
      {
        res.result = colpoPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
    if (resCyto1 == 2) /* LSIL */
    {
      res.result = colpoTEST(COLPOSENSI, location);
      res.nColpo = res.nColpo + 1;
      if (res.result == 0) /* Neg */
      {
        res.result = 2;
      }
      if (res.result == 1) /* CIN 1+ */
      {
        res.result = biopsPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                ColpoBefore);
        res.nBiops = res.nBiops + 1;
      }
    }
    if (resCyto1 == 3) /* HSIL */
    {
      res.result = colpoTEST(COLPOSENSI, location);
      res.nColpo = res.nColpo + 1;
      if (res.result == 0) /* Neg */
      {
        res.result = 4;
      }
      if (res.result == 1) /* Pos */
      {
        res.result = biopsPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                ColpoBefore);
        res.nBiops = res.nBiops + 1;
      }
    }
  }
  return res;
}

struct scr_resultC scrSchema_C(SEXP SCREENSENSI, SEXP SCREENSENSI2, SEXP SCREENSENSI3, SEXP COLPOSENSI,
                               SEXP BIOPSENSI, SEXP HPVTESTSENSI, int location, double probscyto[3],
                                                                                                int cytohpvperiod, int cytohpv_TOCA, int afterSwitch, int hpv_toca,
                                                                                                int cytoType)
{
  int resCyto1=0, resHpvTest=0, ColpoBefore=0;
  struct scr_resultC res={.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  
  if ((cytohpv_TOCA == 1) | (hpv_toca == 1))
  {
    if (hpv_toca == 1)
    {
      res.result = hpvTEST(HPVTESTSENSI, location);
      if (res.result != 0)
      {
        res.result = colpoPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
    if (cytohpv_TOCA==1)
    {
      res.result = cytohpvPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore, probscyto);
      res.nCyto  = res.nCyto + 1;
      res.nHPV   = res.nHPV + 1;
      if (cytoType == 1)
      {
        res.nCombos = res.nCombos + 1;
      }
      if (res.result == 88)
      {
        res.result = 0;
      }
      if ((res.result == 888) | (res.result == 8888))
      {
        res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
  }else{
    if (afterSwitch == 0)
    {
      resCyto1 = cytologyTEST(SCREENSENSI, location, probscyto, 0);
      res.nCyto  = res.nCyto + 1;
      if (resCyto1 == 0)
      {
        res.result = 0;
      }
      if (resCyto1 == 1) /* ASCUS */
      {
        resHpvTest = hpvTEST(HPVTESTSENSI, location);
        res.nHPV   = res.nHPV + 1;
        if (cytoType == 1)
        {
          res.nCombos = res.nCombos + 1;
        }
        if (resHpvTest == 0) /* HPV- or LR */
        {
          res.result = 1;
        }
        if (resHpvTest == 1) /* HPV-HR */
        {
          res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
          res.nColpo = res.nColpo + 1;
          if (res.result == 999)
          {
            res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                    ColpoBefore);
            res.nBiops = res.nBiops + 1;
          }
        }
      }
      if (resCyto1 > 1) /* LSIL+ */
      {
        res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
    
    if (afterSwitch == 1) /* After switch age */
    {
      resHpvTest = hpvTEST(HPVTESTSENSI, location);
      res.nHPV   = res.nHPV + 1;
      if (resHpvTest == 0) /* HPV- or LR */
      {
        res.result = 0;
      }
      if (resHpvTest == 1) /* HPV-HR */
      {
        resCyto1  = cytologyTEST(SCREENSENSI3, location, probscyto, 0);
        res.nCyto = res.nCyto + 1;
        if (cytoType == 1)
        {
          res.nCombos = res.nCombos + 1;
        }
        if (resCyto1 == 0) /* Neg */
        {
          res.result = 2;
        }
        if (resCyto1 > 0) /* ASCUS+ */
        {
          res.result = colpoPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
          res.nColpo = res.nColpo + 1;
          if (res.result == 999)
          {
            res.result = biopsPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                    ColpoBefore);
            res.nBiops = res.nBiops + 1;
          }
        }
      }
    }
  }
  return res;
}

struct scr_resultD scrSchema_D(SEXP SCREENSENSI, SEXP SCREENSENSI2, SEXP SCREENSENSI3, SEXP COLPOSENSI,
                               SEXP BIOPSENSI, SEXP HPVTESTSENSI, int location, double probscyto[3],
                                                                                                int cytohpvperiod, int cytohpv_TOCA, int afterSwitch, int hpv_toca,
                                                                                                int cytoType)
{
  int resCyto1=0, resHpvTest=0, ColpoBefore=0, resHpvGeno=0;
  struct scr_resultD res={.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  
  if ((cytohpv_TOCA == 1) | (hpv_toca == 1))
  {
    if (hpv_toca == 1)
    {
      res.result = hpvTEST(HPVTESTSENSI, location);
      if (res.result != 0)
      {
        res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
    if (cytohpv_TOCA==1)
    {
      res.result = cytohpvPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore, probscyto);
      res.nCyto  = res.nCyto + 1;
      res.nHPV   = res.nHPV + 1;
      if (cytoType == 1)
      {
        res.nCombos = res.nCombos + 1;
      }
      if (res.result == 88)
      {
        res.result = 0;
      }
      if ((res.result == 888) | (res.result == 8888))
      {
        res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
  }else{
    if (afterSwitch == 0)
    {
      resCyto1 = cytologyTEST(SCREENSENSI, location, probscyto, 0);
      res.nCyto  = res.nCyto + 1;
      if (resCyto1 == 0)
      {
        res.result = 0;
      }
      if (resCyto1 == 1) /* ASCUS */
      {
        resHpvTest = hpvTEST(HPVTESTSENSI, location);
        res.nHPV   = res.nHPV + 1;
        if (cytoType == 1)
        {
          res.nCombos = res.nCombos + 1;
        }
        if (resHpvTest == 0) /* HPV- or LR */
        {
          res.result = 1;
        }
        if (resHpvTest == 1) /* HPV-HR */
        {
          res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
          res.nColpo = res.nColpo + 1;
          if (res.result == 999)
          {
            res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                    ColpoBefore);
            res.nBiops = res.nBiops + 1;
          }
        }
      }
      if (resCyto1 > 1) /* LSIL+ */
      {
        res.result = colpoPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
        res.nColpo = res.nColpo + 1;
        if (res.result == 999)
        {
          res.result = biopsPROTO(SCREENSENSI2, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                  ColpoBefore);
          res.nBiops = res.nBiops + 1;
        }
      }
    }
    
    if (afterSwitch == 1) /* After switch age */
    {
      resHpvTest = hpvTEST(HPVTESTSENSI, location);
      res.nHPV   = res.nHPV + 1;
      if (resHpvTest == 0) /* HPV- or LR */
      {
        res.result = 0;
      }
      if (resHpvTest == 1) /* HPV-HR */
      {
        resHpvGeno = rbinom(1, 0.3); /* 16/18 70% */
        if (resHpvGeno == 0) /* HPV 16/18 */
        {
          res.result = colpoPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
          res.nColpo = res.nColpo + 1;
          if (res.result == 999)
          {
            res.result = biopsPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                    ColpoBefore);
            res.nBiops = res.nBiops + 1;
          }
        }
        if (resHpvGeno == 1) /* HPV other HR */
        {
          resCyto1 = cytologyTEST(SCREENSENSI3, location, probscyto, 0);
          res.nCyto = res.nCyto + 1;
          if (cytoType == 1)
          {
            res.nCombos = res.nCombos + 1;
          }
          if (resCyto1 == 0) /* Neg */
          {
            res.result = 2;
          }
          if (resCyto1 > 0) /* ASCUS+ */
          {
            res.result = colpoPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, location, ColpoBefore);
            res.nColpo = res.nColpo + 1;
            if (res.result == 999)
            {
              res.result = biopsPROTO(SCREENSENSI3, HPVTESTSENSI, COLPOSENSI, BIOPSENSI, location,
                                      ColpoBefore);
              res.nBiops = res.nBiops + 1;
            }
          }
        }
      }
    }
  }
  return res;
}

SEXP Cmicrosim(SEXP TRANSITION, SEXP NSTATES, SEXP ABS_STATES, SEXP NABS_STATES, SEXP SYMPT_STATES, SEXP NSYMPTSTATES, SEXP PROB_SYMPT,
               SEXP SIZE, SEXP P_MEN, SEXP STEPS, SEXP MIN_AGE, SEXP DAGE, SEXP UTILITYCOEFS, SEXP COSTCOEFS_MD,
               SEXP COSTCOEFS_NMD, SEXP COSTCOEFS_I, SEXP DISC, SEXP VACC, SEXP VACC_AGE, SEXP NDOSES, SEXP VACC_COV,
               SEXP VACC_EFF, SEXP VACCPRICE_MD, SEXP VACCPRICE_NMD, SEXP VACCPRICE_I, SEXP NVACCPERIODS,
               SEXP SCREENING, SEXP SCREENTYPE, SEXP SCRSCHEMA, SEXP SCREENPERIOD, SEXP CYTOTYPE,
               SEXP SCREENPRICE_MD, SEXP SCREENPRICE_NMD, SEXP SCREENPRICE_I, SEXP COLPOPRICE_MD, SEXP COLPOPRICE_NMD,
               SEXP COLPOPRICE_I, SEXP BIOPSPRICE_MD, SEXP BIOPSPRICE_NMD, SEXP BIOPSPRICE_I, SEXP SCREENCOVERAGE,
               SEXP SCREENSENSI, SEXP SCREENSENSI2, SEXP SCREENSENSI3, SEXP COLPOSENSI, SEXP BIOPSENSI,
               SEXP CYTOHPVPRICE_MD, SEXP CYTOHPVPRICE_NMD, SEXP CYTOHPVPRICE_I, SEXP HPVTESTSENSI,
			         SEXP HPVTESTPRICE_MD, SEXP HPVTESTPRICE_NMD, SEXP HPVTESTPRICE_I, SEXP TREATPROBS,
			         SEXP NANNUALVISITS, SEXP NANNUALVISITSLSIL, SEXP NANNUALVISITSHSIL, SEXP CYTOHPVPERIOD,
			         SEXP CYTOHPVPOSTCOLPO, SEXP CYTOHPVPOSTBIOP, SEXP CYTOLSILPERIOD, SEXP CYTOHSILPERIOD, SEXP SWITCHAGE,
			         SEXP C_PERIOD, SEXP HPVPERIOD)
{
  int i=0, j=0, k=0, l=0, c=0, ind=0, *dage=INTEGER(DAGE), *size=INTEGER(SIZE),
    *steps=INTEGER(STEPS), *nstates=INTEGER(NSTATES), *min_age=INTEGER(MIN_AGE),
    *nsymptstates=INTEGER(NSYMPTSTATES), ageG, row, location=0, maximum, symptstateB=0, ageStep,
    *vacc=INTEGER(VACC), *nVaccPeriods=INTEGER(NVACCPERIODS), *screen=INTEGER(SCREENING),
    *scrSchema=INTEGER(SCRSCHEMA), *scrPeriod=INTEGER(SCREENPERIOD), screenedW=0, *cytoType=INTEGER(CYTOTYPE),
    curedW=0, change_state=0, vaccStep=0, vaccinatedW=0, indVacc=0, ageVacc=0, treatedW=0, ind2=0, ind3=0, ind4=0,
    ind5=0, ind6=0, ind7=0, ind8=0, ind9=0, pos=0,
    *nAnnualVisits=INTEGER(NANNUALVISITS), *nAnnualVisitsLSIL=INTEGER(NANNUALVISITSLSIL), *nAnnualVisitsHSIL=INTEGER(NANNUALVISITSHSIL),
    *cytohpvperiod=INTEGER(CYTOHPVPERIOD), *cytohpvpostcolpo=INTEGER(CYTOHPVPOSTCOLPO),
    *cytoLSILperiod=INTEGER(CYTOLSILPERIOD), *cytoHSILperiod=INTEGER(CYTOHSILPERIOD), pathPROTO=0,
	  cytohpv_TOCA=0, cyto_TOCA=0, *switchAge=INTEGER(SWITCHAGE), afterSwitch=0, hpv_toca=0,
          *hpvPeriod=INTEGER(HPVPERIOD), newCytos[*steps], newColpo[*steps], newHPVTests[*steps],
         newBiops[*steps], CytoW=0, HPVTestW=0, ColpoW=0, BiopW=0, SymptW=0, TreatedW2=0, nCytoStep=0,
         nColpoStep=0, nHPVTestStep=0, nBiopsyStep=0;
  double *p_men=REAL(P_MEN), *disc=REAL(DISC), md_cost[*size], *costCoefs_md, nmd_cost[*size], i_cost[*size],
         *costCoefs_nmd, *costCoefs_i, utilities[*size], *utility, md_cost_disc[*size], newHPV[*steps], newCC[*steps],
         newCIN1[*steps], newCIN2[*steps], newCIN3[*steps], newCCMDeath[*steps], newOtherDeath[*steps],
         nmd_cost_disc[*size], i_cost_disc[*size], utilities_disc[*size], le[*size], le_disc[*size],
         *vaccPrice_md=REAL(VACCPRICE_MD), *vaccPrice_nmd=REAL(VACCPRICE_NMD),
         *vaccPrice_i=REAL(VACCPRICE_I), *screenPrice_md=REAL(SCREENPRICE_MD), *screenPrice_nmd=REAL(SCREENPRICE_NMD),
         *screenPrice_i=REAL(SCREENPRICE_I), *colpoPrice_md=REAL(COLPOPRICE_MD), *colpoPrice_nmd=REAL(COLPOPRICE_NMD),
         *colpoPrice_i=REAL(COLPOPRICE_I), *hpvTestPrice_md=REAL(HPVTESTPRICE_MD), *hpvTestPrice_nmd=REAL(HPVTESTPRICE_NMD),
         *hpvTestPrice_i=REAL(HPVTESTPRICE_I), *biopsPrice_md=REAL(BIOPSPRICE_MD), *biopsPrice_nmd=REAL(BIOPSPRICE_NMD),
         *biopsPrice_i=REAL(BIOPSPRICE_I), *cytoHpvPrice_md=REAL(CYTOHPVPRICE_MD), *cytoHpvPrice_nmd=REAL(CYTOHPVPRICE_NMD),
         *cytoHpvPrice_i=REAL(CYTOHPVPRICE_I), probscyto[3]={0.73, 0.22, 0.05}; /* Int J Biomed Sci. 2013 Dec; 9(4): 205–210 */
  struct scr_resultA global_scrA = {.nCyto=0, .nColpo=0, .result=0};
  struct scr_resultB global_scrB = {.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  struct scr_resultC global_scrC = {.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  struct scr_resultD global_scrD = {.nCyto=0, .nHPV=0, .nColpo=0, .nBiops=0, .nCombos=0, .result=0};
  SEXP RES;
  SEXP RES2;
  SEXP r;
  SEXP trans_v;
  SEXP res_scr;
  PROTECT(RES     = allocMatrix(REALSXP, *size+11, *steps+18));
  PROTECT(RES2    = allocVector(INTSXP, *nstates));
  PROTECT(r       = allocMatrix(REALSXP, *nstates, *nstates));
  PROTECT(trans_v = allocVector(REALSXP, *nstates));
  PROTECT(res_scr = allocVector(INTSXP, 1));

  GetRNGstate();

  for (i=0; i<*steps; i++)
  {
    newHPV[i]        = 0.0;
    newCC[i]         = 0.0;
    newCIN1[i]       = 0.0;
    newCIN2[i]       = 0.0;
    newCIN3[i]       = 0.0;
    newCCMDeath[i]   = 0.0;
    newOtherDeath[i] = 0.0;
    newCytos[i]      = 0;
    newColpo[i]      = 0;
    newHPVTests[i]   = 0;
    newBiops[i]      = 0;
  }
  for (i=0; i<*size; i++)
  {
    vaccinatedW       = 0;
    ageG              = 0;
    REAL(RES)[i]      = 0.0;
    md_cost[i]        = 0.0;
    nmd_cost[i]       = 0.0;
    i_cost[i]         = 0.0;
    md_cost_disc[i]   = 0.0;
    nmd_cost_disc[i]  = 0.0;
    i_cost_disc[i]    = 0.0;
    utilities[i]      = 0.5; /*1.0 */
    utilities_disc[i] = 0.5; /*1.0 */
    le[i]             = 0.5; /*1.0 */
    le_disc[i]        = 0.5; /*1.0 */
    ind2              = 0;
    CytoW             = 0;
    HPVTestW          = 0;
    ColpoW            = 0;
    BiopW             = 0;
    SymptW            = 0;
    TreatedW2         = 0;
    cyto_TOCA         = 0;
    cytohpv_TOCA      = 0;
    hpv_toca          = 0;
    ind3              = 0;
    ind4              = 0;
    ind5              = 0;
    ind6              = 0;
    ind7              = 0;
    ind8              = 0;
    ind9              = 0;
    for (j=1; j<*steps; j++)
    {
      change_state = 1;
      /*ageStep      = *min_age + j;*/ /* Step is 1 year */
      ageStep      = (int)(*min_age + j/2); /* Step is 1 semester */
      screenedW    = 0;
      treatedW     = 0;
      curedW       = 0;
      pos          = 0;
      global_scrA.nCyto=0;
      global_scrA.nColpo=0;
      global_scrB.nCyto=0;
      global_scrB.nHPV=0;
      global_scrB.nColpo=0;
      global_scrB.nBiops=0;
      global_scrC.nCyto=0;
      global_scrC.nHPV=0;
      global_scrC.nColpo=0;
      global_scrC.nBiops=0;
      global_scrD.nCyto=0;
      global_scrD.nHPV=0;
      global_scrD.nColpo=0;
      global_scrD.nBiops=0;
      nCytoStep=0;
      nColpoStep=0;
      nHPVTestStep=0;
      nBiopsyStep=0;
      if ((*screen == 1) & ((*scrSchema == 3) | (*scrSchema == 4)))
      {
        if (ageStep <= *switchAge)
        {
          afterSwitch = 0;
          *scrPeriod  = *INTEGER(VECTOR_ELT(C_PERIOD, 0));
        }
        if (ageStep > *switchAge)
        {
          afterSwitch = 1;
          *scrPeriod  = *INTEGER(VECTOR_ELT(C_PERIOD, 1));
        }
      }
      if (((j % 2 == 0)) & (((j/2) % *dage) == 0))
      {
        ageG = ageG+1;
      }
      if (REAL(RES)[i+(*size+11)*(j-1)] >= *INTEGER(VECTOR_ELT(ABS_STATES, 0))) /* absorbing states: should be the last ones in the transition matrix */
      {
        REAL(RES)[i+(*size+11)*j] = REAL(RES)[i+(*size+11)*(j-1)];
       }else{
        row = REAL(RES)[i+(j-1)*(*size+11)]; /* Previous health state */
        for (k=0; k<*nstates; k++)
        {
          for (l=0; l<*nstates; l++)
          {
            REAL(r)[k+*nstates*l] = REAL(VECTOR_ELT(TRANSITION, l))[k+(*nstates)*ageG];
            REAL(trans_v)[l]      = REAL(r)[row+(*nstates)*l];
          }
        }
        /* vaccination module */
        vaccStep = 0;
        indVacc  = 0;
        if ((*vacc == 1) & (vaccinatedW == 0))
        {
          for (ind=0; ind<*nVaccPeriods; ind++)
          {
            if (*REAL(VECTOR_ELT(VACC_AGE, ind)) == ageStep)
            {
              vaccStep = 1;
              indVacc = ind;
            }
          }
          if ((*vacc == 1) & (vaccStep == 1))
          {
            vaccinatedW = rbinom(1, *REAL(VECTOR_ELT(VACC_COV, indVacc)));
            if (vaccinatedW == 1)
            {
              ageVacc          = indVacc;
              md_cost[i]       = md_cost[i]       + *vaccPrice_md*(*REAL(VECTOR_ELT(NDOSES, indVacc))); /* vaccination costs */
              nmd_cost[i]      = nmd_cost[i]      + *vaccPrice_nmd*(*REAL(VECTOR_ELT(NDOSES, indVacc)));
              i_cost[i]        = i_cost[i]        + *vaccPrice_i*(*REAL(VECTOR_ELT(NDOSES, indVacc)));
              md_cost_disc[i]  = md_cost_disc[i]  + *vaccPrice_md*(*REAL(VECTOR_ELT(NDOSES, indVacc)))/pow((1+*disc/200), (j-1)); /* vaccination costs */
              nmd_cost_disc[i] = nmd_cost_disc[i] + *vaccPrice_nmd*(*REAL(VECTOR_ELT(NDOSES, indVacc)))/pow((1+*disc/200), (j-1));
              i_cost_disc[i]   = i_cost_disc[i]   + *vaccPrice_i*(*REAL(VECTOR_ELT(NDOSES, indVacc)))/pow((1+*disc/200), (j-1));
            }
          }
        }
        if ((vaccinatedW == 1) & (row == 0)) /* modify HPV infection probabilities if vaccinated */
        {
          REAL(trans_v)[1] = (1.0-*REAL(VECTOR_ELT(VACC_EFF, ageVacc)))*REAL(trans_v)[1];
          REAL(trans_v)[0] = 1.0-REAL(trans_v)[1]-REAL(trans_v)[2]-REAL(trans_v)[3]-
            REAL(trans_v)[4]-REAL(trans_v)[5]-REAL(trans_v)[6]-REAL(trans_v)[7]-
            REAL(trans_v)[8]-REAL(trans_v)[9]-REAL(trans_v)[10]-REAL(trans_v)[11];
        }
        /*end of vaccination module */
        rmultinom(1, REAL(trans_v), *nstates, INTEGER(RES2)); /* set current state */

        maximum = INTEGER(RES2)[0];
        location = 0;
        for (c=1; c<*nstates; c++)
        {
          if (INTEGER(RES2)[c] >= maximum)
          {
            maximum  = INTEGER(RES2)[c];
            location = c;
          }
       }
       REAL(RES)[i+(*size+11)*j] = location;
       utility       = REAL(VECTOR_ELT(UTILITYCOEFS, location));
       costCoefs_md  = REAL(VECTOR_ELT(COSTCOEFS_MD, location));
       costCoefs_nmd = REAL(VECTOR_ELT(COSTCOEFS_NMD, location));
       costCoefs_i   = REAL(VECTOR_ELT(COSTCOEFS_I, location));

       if ((location != 10) & (location != 11))
       {
         utilities[i]      = utilities[i]      + *utility/2;
         utilities_disc[i] = utilities_disc[i] + (*utility/2)/pow((1+*disc/200), (j-1));
         le[i]             = le[i] + 0.5;
         le_disc[i]        = le_disc[i] + 0.5/pow((1+*disc/100), *min_age+0.5);
         /* ORGANIZED SCREENING MODULE */
         if (*screen == 1)
         {
           if (2*j % *scrPeriod == 0)
           {
             screenedW = rbinom(1, *REAL(VECTOR_ELT(SCREENCOVERAGE, ageG)));
           }
           if ((*scrSchema==1) & (ind2>0))
           {
             screenedW = 1;
           }
           if (*scrSchema==2)
           {
           	 if (((ind4 == *cytoLSILperiod) & (ind4>0) & (ind8 < *nAnnualVisitsLSIL)) | ((ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 0))) & (ind5>0)) |
              ((ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) & (ind5>0)) | ((ind6 == *cytohpvpostcolpo) & (ind6>0)) |
              ((ind7 == *cytoHSILperiod) & (ind7>0) & ((ind9 < *nAnnualVisitsHSIL) & (ind9>0))))
             {
               screenedW = 1;
             }
           }
           if ((*scrSchema==3) | (*scrSchema==4))
           {
           	 if (((ind3 == *cytohpvperiod) & (ind3>0)) | ((ind4 == *hpvPeriod) & (ind4>0)) | ((ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 0))) & (ind5>0)) |
                 ((ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) & (ind5>0)) | ((ind6 == *cytohpvpostcolpo) & (ind6>0)))
             {
               screenedW = 1;
             }
           }
           if ((*scrSchema == 1) & (screenedW == 1))
         	 {
         	   if ((ind2 >= *nAnnualVisits) & (ind2 != 999))
         	   {
         	 	   ind2 = 0;
         	   }
         	   if (ind2 > 0)
             {
           	   if (ind2 == 999)
           	   {
           	 	    ind2 = 1;
           	   }else{
           	 	    ind2 = ind2 + 1;
           	   }
             }
         	 }
         	 if ((*scrSchema == 2) & (screenedW == 1))
         	 {
         	   if (((ind4 > *cytoLSILperiod) & (ind4 != 999)) | (ind8 >= *nAnnualVisitsLSIL))
         	   {
         	 	   ind4 = 0;
               ind8 = 0;
         	   }
         	   if ((ind5 > *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) & (ind5 != 999))
             {
               ind5 = 0;
             }
             if ((ind6 > *cytohpvpostcolpo) & (ind6 != 999))
             {
               ind6 = 0;
             }
         	   if (((ind7 > *cytoHSILperiod) & (ind7 != 999)) | (ind9 >= *nAnnualVisitsHSIL))
         	   {
         	 	   ind7 = 0;
               ind9 = 0;
         	   }
         	   if (ind4 > 0)
             {
           	   if (ind4 == 999)
           	   {
           	     ind4 = 1;
								 ind8 = 1;
           	   }else{
           	 	   ind4 = ind4 + 1;
                 ind8 = ind8 + 1;
               }
             }
             if (ind5 > 0)
             {
           	   if (ind5 == 999)
           	   {
                 ind5 = 1;
               }else{
               	 ind5 = ind5 + 1;
               }
             }
             if (ind6 > 0)
             {
           	   if (ind6 == 999)
           	   {
                 ind6 = 1;
               }else{
             	   ind6 = ind6 + 1;
               }
             }
             if (ind7 > 0)
             {
           	   if (ind7 == 999)
           	   {
                 ind7 = 1;
								 ind9 = 1;
               }else{
             	   ind7 = ind7 + 1;
                 ind9 = ind9 + 1;
               }
             }
         	 }
         	 if (((*scrSchema == 3) | (*scrSchema == 4)) & (screenedW == 1))
         	 {
         	   if ((ind3 > *cytohpvperiod) & (ind3 != 999))
             {
               ind3 = 0;
             }
             if ((ind4 > *hpvPeriod) & (ind4 != 999))
             {
               ind4 = 0;
             }
             if ((ind5 > *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) & (ind5 != 999))
             {
               ind5 = 0;
             }
             if ((ind6 > *cytohpvpostcolpo) & (ind6 != 999))
             {
               ind6 = 0;
             }
         	   if (ind3 > 0)
             {
           	   if (ind3 == 999)
           	   {
                 ind3 = 1;
               }else{
                 ind3 = ind3 + 1;
               }
             }
         	   if (ind4 > 0)
             {
           	   if (ind4 == 999)
           	   {
           	     ind4 = 1;
           	   }else{
           	 	   ind4 = ind4 + 1;
               }
             }
             if (ind5 > 0)
             {
           	   if (ind5 == 999)
           	   {
                 ind5 = 1;
               }else{
               	 ind5 = ind5 + 1;
               }
             }
             if (ind6 > 0)
             {
           	   if (ind6 == 999)
           	   {
                 ind6 = 1;
               }else{
             	   ind6 = ind6 + 1;
               }
             }
           }
         }

         if ((screenedW == 1) & (*scrSchema == 1)) /* SCHEMA A: Cytology alone with repeat cytology for triage */
         {
           if (ind2>0)
           {
             pos = 1;
           }else{
           	 pos = 0;
           }
           global_scrA = scrSchema_A(SCREENSENSI, SCREENSENSI2, COLPOSENSI, location, pos, probscyto);
           nCytoStep   = nCytoStep   + global_scrA.nCyto;
           newCytos[j] = newCytos[j] + nCytoStep;
           nColpoStep  = nColpoStep  + global_scrA.nColpo;
           newColpo[j] = newColpo[j] + nColpoStep;
           if (global_scrA.result == 2)
           {
             ind2 = 999;
           }
           if (global_scrA.result == 3)
           {
             treatedW = 1;
             ind2     = 0;
             pos      = 0;
           }
         }

         if ((screenedW == 1) & (*scrSchema == 2)) /* SCHEMA B: Cytology with HPV triage */
         {
         	 if ((ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 0))) | (ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) |
               (ind6 == *cytohpvpostcolpo))
           {
             cytohpv_TOCA = 1;
           }else{
           	 cytohpv_TOCA = 0;
           }
           if (((ind4 == *cytoLSILperiod) & (ind8 != 0)) | ((ind7 == *cytoHSILperiod) & (ind9 != 0)))
           {
             cyto_TOCA = 1;
           }else{
           	 cyto_TOCA = 0;
           }
           global_scrB = scrSchema_B(SCREENSENSI, SCREENSENSI2, SCREENSENSI3, COLPOSENSI, BIOPSENSI, HPVTESTSENSI,
             location, probscyto, pathPROTO, cytohpv_TOCA, cyto_TOCA, *cytoType);
           nCytoStep      = nCytoStep      + global_scrB.nCyto;
           newCytos[j]    = newCytos[j]    + nCytoStep;
           nColpoStep     = nColpoStep     + global_scrB.nColpo;
           newColpo[j]    = newColpo[j]    + nColpoStep;
           nHPVTestStep   = nHPVTestStep   + global_scrB.nHPV;
           newHPVTests[j] = newHPVTests[j] + nHPVTestStep;
           nBiopsyStep    = nBiopsyStep    + global_scrB.nBiops;
           newBiops[j]    = newBiops[j]    + nBiopsyStep;
           if (global_scrB.result == 2)
           {
             ind4 = 999;
           }
           if (global_scrB.result == 77)
           {
             ind5 = 999;
           }
           if (global_scrB.result == 99)
           {
             ind6 = 999;
           }
           if (global_scrB.result == 4)
           {
             ind7 = 999;
           }
           if (global_scrB.result == 777)
           {
             treatedW = 1;
             ind4     = 0;
             ind5     = 0;
             ind6     = 0;
             ind7     = 0;
             ind8     = 0;
             ind9     = 0;
           }
         }

         if ((screenedW == 1) & (*scrSchema == 3))  /* SCHEMA C: HPV with cytology triage */
         {
           if ((ind3 == *cytohpvperiod) | (ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 0))) | (ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) |
               (ind6 == *cytohpvpostcolpo))
           {
             cytohpv_TOCA = 1;
           }else{
           	 cytohpv_TOCA = 0;
           }
           if (ind4 == *hpvPeriod)
           {
             hpv_toca = 1;
           }else{
           	 hpv_toca = 0;
           }
           global_scrC = scrSchema_C(SCREENSENSI, SCREENSENSI2, SCREENSENSI3, COLPOSENSI,
                                     BIOPSENSI, HPVTESTSENSI, location, probscyto, *cytohpvperiod, cytohpv_TOCA, afterSwitch, hpv_toca, *cytoType);
           nCytoStep   = nCytoStep + global_scrC.nCyto;
           newCytos[j] = newCytos[j] + nCytoStep;
           nColpoStep  = nColpoStep + global_scrC.nColpo;
           newColpo[j] = newColpo[j] + nColpoStep;
           nHPVTestStep = nHPVTestStep + global_scrC.nHPV;
           newHPVTests[j] = newHPVTests[j] + nHPVTestStep;
           nBiopsyStep = nBiopsyStep + global_scrC.nBiops;
           newBiops[j] = newBiops[j] + nBiopsyStep;
           if (global_scrC.result == 1)
           {
             ind3 = 999;
           }
           if (global_scrC.result == 2)
           {
             ind4 = 999;
           }
           if (global_scrC.result == 77)
           {
             ind5 = 999;
           }
           if (global_scrC.result == 99)
           {
             ind6 = 999;
           }
           if (global_scrC.result == 777)
           {
             treatedW = 1;
             ind3     = 0;
             ind4     = 0;
             ind5     = 0;
             ind6     = 0;
           }
         }
         if ((screenedW == 1) & (*scrSchema == 4))  /* SCHEMA D: HPV genotyping with cytology triage */
         {
         	 if ((ind3 == *cytohpvperiod) | (ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 0))) | 
         	 	   (ind5 == *INTEGER(VECTOR_ELT(CYTOHPVPOSTBIOP, 1))) |
               (ind6 == *cytohpvpostcolpo))
           {
             cytohpv_TOCA = 1;
           }else{
           	 cytohpv_TOCA = 0;
           }
           if (ind4 == *hpvPeriod)
           {
             hpv_toca = 1;
           }else{
           	 hpv_toca = 0;
           }
           global_scrD = scrSchema_D(SCREENSENSI, SCREENSENSI2, SCREENSENSI3, COLPOSENSI,
                                     BIOPSENSI, HPVTESTSENSI, location, probscyto, *cytohpvperiod, cytohpv_TOCA, afterSwitch, hpv_toca, *cytoType);
           nCytoStep   = nCytoStep + global_scrD.nCyto;
           newCytos[j] = newCytos[j] + nCytoStep;
           nColpoStep  = nColpoStep + global_scrD.nColpo;
           newColpo[j] = newColpo[j] + nColpoStep;
           nHPVTestStep = nHPVTestStep + global_scrD.nHPV;
           newHPVTests[j] = newHPVTests[j] + nHPVTestStep;
           nBiopsyStep = nBiopsyStep + global_scrD.nBiops;
           newBiops[j] = newBiops[j] + nBiopsyStep;
           if (global_scrD.result == 1)
           {
             ind3 = 999;
           }
           if (global_scrD.result == 2)
           {
             ind4 = 999;
           }
           if (global_scrD.result == 77)
           {
             ind5 = 999;
           }
           if (global_scrD.result == 99)
           {
             ind6 = 999;
           }
           if (global_scrD.result == 777)
           {
             treatedW = 1;
             ind3     = 0;
             ind4     = 0;
             ind5     = 0;
             ind6     = 0;
           }
         }

         if ((*screen == 1) & (*scrSchema == 1))
         {
           md_cost[i]   = md_cost[i]  + global_scrA.nCyto*(*screenPrice_md)  + global_scrA.nColpo*(*colpoPrice_md);
           nmd_cost[i]  = nmd_cost[i] + global_scrA.nCyto*(*screenPrice_nmd) + global_scrA.nColpo*(*colpoPrice_nmd);
           i_cost[i]    = i_cost[i]   + global_scrA.nCyto*(*screenPrice_i)   + global_scrA.nColpo*(*colpoPrice_i);

           md_cost_disc[i]   = md_cost_disc[i]  + (global_scrA.nCyto*(*screenPrice_md)  +
             global_scrA.nColpo*(*colpoPrice_md))/pow((1+*disc/200), (j-1));
           nmd_cost_disc[i]  = nmd_cost_disc[i] + (global_scrA.nCyto*(*screenPrice_nmd) +
             global_scrA.nColpo*(*colpoPrice_nmd))/pow((1+*disc/200), (j-1));
           i_cost_disc[i]    = i_cost_disc[i]   + (global_scrA.nCyto*(*screenPrice_i)   +
             global_scrA.nColpo*(*colpoPrice_i))/pow((1+*disc/200), (j-1));
         }
         if ((*screen == 1) & (*scrSchema == 2))
         {
           md_cost[i]   = md_cost[i]  + global_scrB.nCyto*(*screenPrice_md)  + global_scrB.nHPV*(*hpvTestPrice_md) +
             global_scrB.nColpo*(*colpoPrice_md) + global_scrB.nBiops*(*biopsPrice_md) +
             global_scrB.nCombos*(*cytoHpvPrice_md) - global_scrB.nCombos*(*screenPrice_md) -
             global_scrB.nCombos*(*hpvTestPrice_md);
           nmd_cost[i]  = nmd_cost[i] + global_scrB.nCyto*(*screenPrice_nmd) + global_scrB.nHPV*(*hpvTestPrice_nmd) +
             global_scrB.nColpo*(*colpoPrice_nmd) + global_scrB.nBiops*(*biopsPrice_nmd) +
             global_scrB.nCombos*(*cytoHpvPrice_nmd) - global_scrB.nCombos*(*screenPrice_nmd) -
             global_scrB.nCombos*(*hpvTestPrice_nmd);
           i_cost[i]    = i_cost[i]   + global_scrB.nCyto*(*screenPrice_i)   + global_scrB.nHPV*(*hpvTestPrice_i) +
             global_scrB.nColpo*(*colpoPrice_i) + global_scrB.nBiops*(*biopsPrice_i) +
             global_scrB.nCombos*(*cytoHpvPrice_i) - global_scrB.nCombos*(*screenPrice_i) -
             global_scrB.nCombos*(*hpvTestPrice_i);

           md_cost_disc[i]   = md_cost_disc[i]  + (global_scrB.nCyto*(*screenPrice_md) +
             global_scrB.nHPV*(*hpvTestPrice_md) + global_scrB.nColpo*(*colpoPrice_md) +
             global_scrB.nBiops*(*biopsPrice_md) + global_scrB.nCombos*(*cytoHpvPrice_md) -
             global_scrB.nCombos*(*screenPrice_md) - global_scrB.nCombos*(*hpvTestPrice_md))/pow((1+*disc/200), (j-1));
           nmd_cost_disc[i]  = nmd_cost_disc[i] + (global_scrB.nCyto*(*screenPrice_nmd) +
             global_scrB.nHPV*(*hpvTestPrice_nmd) + global_scrB.nColpo*(*colpoPrice_nmd) +
             global_scrB.nBiops*(*biopsPrice_nmd) + global_scrB.nCombos*(*cytoHpvPrice_nmd) -
             global_scrB.nCombos*(*screenPrice_nmd) - global_scrB.nCombos*(*hpvTestPrice_nmd))/pow((1+*disc/200), (j-1));
           i_cost_disc[i]    = i_cost_disc[i]   + (global_scrB.nCyto*(*screenPrice_i) +
             global_scrB.nHPV*(*hpvTestPrice_i) + global_scrB.nColpo*(*colpoPrice_i) +
             global_scrB.nBiops*(*biopsPrice_i) + global_scrB.nCombos*(*cytoHpvPrice_i) -
             global_scrB.nCombos*(*screenPrice_i) - global_scrB.nCombos*(*hpvTestPrice_i))/pow((1+*disc/200), (j-1));
         }
         if ((*screen == 1) & (*scrSchema == 3))
         {
           md_cost[i]   = md_cost[i]  + global_scrC.nCyto*(*screenPrice_md)  + global_scrC.nHPV*(*hpvTestPrice_md) +
             global_scrC.nColpo*(*colpoPrice_md) + global_scrC.nBiops*(*biopsPrice_md) +
             global_scrC.nCombos*(*cytoHpvPrice_md) - global_scrC.nCombos*(*screenPrice_md) -
             global_scrC.nCombos*(*hpvTestPrice_md);
           nmd_cost[i]  = nmd_cost[i] + global_scrC.nCyto*(*screenPrice_nmd) + global_scrC.nHPV*(*hpvTestPrice_nmd) +
             global_scrC.nColpo*(*colpoPrice_nmd) + global_scrC.nBiops*(*biopsPrice_nmd) +
             global_scrC.nCombos*(*cytoHpvPrice_nmd) - global_scrC.nCombos*(*screenPrice_nmd) -
             global_scrC.nCombos*(*hpvTestPrice_nmd);
           i_cost[i]    = i_cost[i]   + global_scrC.nCyto*(*screenPrice_i)   + global_scrC.nHPV*(*hpvTestPrice_i) +
             global_scrC.nColpo*(*colpoPrice_i) + global_scrC.nBiops*(*biopsPrice_i) +
             global_scrC.nCombos*(*cytoHpvPrice_i) - global_scrC.nCombos*(*screenPrice_i) -
             global_scrC.nCombos*(*hpvTestPrice_i);

           md_cost_disc[i]   = md_cost_disc[i]  + (global_scrC.nCyto*(*screenPrice_md) +
             global_scrC.nHPV*(*hpvTestPrice_md) + global_scrC.nColpo*(*colpoPrice_md) +
             global_scrC.nBiops*(*biopsPrice_md) + global_scrC.nCombos*(*cytoHpvPrice_md) -
             global_scrC.nCombos*(*screenPrice_md) - global_scrC.nCombos*(*hpvTestPrice_md))/pow((1+*disc/200), (j-1));
           nmd_cost_disc[i]  = nmd_cost_disc[i] + (global_scrC.nCyto*(*screenPrice_nmd) +
             global_scrC.nHPV*(*hpvTestPrice_nmd) + global_scrC.nColpo*(*colpoPrice_nmd) +
             global_scrC.nBiops*(*biopsPrice_nmd) + global_scrC.nCombos*(*cytoHpvPrice_nmd) -
             global_scrC.nCombos*(*screenPrice_nmd) - global_scrC.nCombos*(*hpvTestPrice_nmd))/pow((1+*disc/200), (j-1));
           i_cost_disc[i]    = i_cost_disc[i]   + (global_scrC.nCyto*(*screenPrice_i) +
             global_scrC.nHPV*(*hpvTestPrice_i) + global_scrC.nColpo*(*colpoPrice_i) +
             global_scrC.nBiops*(*biopsPrice_i) + global_scrC.nCombos*(*cytoHpvPrice_i) -
             global_scrC.nCombos*(*screenPrice_i) - global_scrC.nCombos*(*hpvTestPrice_i))/pow((1+*disc/200), (j-1));
         }
         if ((*screen == 1) & (*scrSchema == 4))
         {
           md_cost[i]   = md_cost[i]  + global_scrD.nCyto*(*screenPrice_md)  + global_scrD.nHPV*(*hpvTestPrice_md) +
             global_scrD.nColpo*(*colpoPrice_md) + global_scrD.nBiops*(*biopsPrice_md) +
             global_scrD.nCombos*(*cytoHpvPrice_md) - global_scrD.nCombos*(*screenPrice_md) -
             global_scrD.nCombos*(*hpvTestPrice_md);
           nmd_cost[i]  = nmd_cost[i] + global_scrD.nCyto*(*screenPrice_nmd) + global_scrD.nHPV*(*hpvTestPrice_nmd) +
             global_scrD.nColpo*(*colpoPrice_nmd) + global_scrD.nBiops*(*biopsPrice_nmd) +
             global_scrD.nCombos*(*cytoHpvPrice_nmd) - global_scrD.nCombos*(*screenPrice_nmd) -
             global_scrD.nCombos*(*hpvTestPrice_nmd);
           i_cost[i]    = i_cost[i]   + global_scrD.nCyto*(*screenPrice_i)   + global_scrD.nHPV*(*hpvTestPrice_i) +
             global_scrD.nColpo*(*colpoPrice_i) + global_scrD.nBiops*(*biopsPrice_i) +
             global_scrD.nCombos*(*cytoHpvPrice_i) - global_scrD.nCombos*(*screenPrice_i) -
             global_scrD.nCombos*(*hpvTestPrice_i);

           md_cost_disc[i]   = md_cost_disc[i]  + (global_scrD.nCyto*(*screenPrice_md) +
             global_scrD.nHPV*(*hpvTestPrice_md) + global_scrD.nColpo*(*colpoPrice_md) +
             global_scrD.nBiops*(*biopsPrice_md) + global_scrD.nCombos*(*cytoHpvPrice_md) -
             global_scrD.nCombos*(*screenPrice_md) - global_scrD.nCombos*(*hpvTestPrice_md))/pow((1+*disc/200), (j-1));
           nmd_cost_disc[i]  = nmd_cost_disc[i] + (global_scrD.nCyto*(*screenPrice_nmd) +
             global_scrD.nHPV*(*hpvTestPrice_nmd) + global_scrD.nColpo*(*colpoPrice_nmd) +
             global_scrD.nBiops*(*biopsPrice_nmd) + global_scrD.nCombos*(*cytoHpvPrice_nmd) -
             global_scrD.nCombos*(*screenPrice_nmd) - global_scrD.nCombos*(*hpvTestPrice_nmd))/pow((1+*disc/200), (j-1));
           i_cost_disc[i]    = i_cost_disc[i]   + (global_scrD.nCyto*(*screenPrice_i) +
             global_scrD.nHPV*(*hpvTestPrice_i) + global_scrD.nColpo*(*colpoPrice_i) +
             global_scrD.nBiops*(*biopsPrice_i) + global_scrD.nCombos*(*cytoHpvPrice_i) -
             global_scrD.nCombos*(*screenPrice_i) - global_scrD.nCombos*(*hpvTestPrice_i))/pow((1+*disc/200), (j-1));
         }

         symptstateB = 0;
         for (c=0; c<*nsymptstates; c++)
         {
           if (location == *REAL(VECTOR_ELT(SYMPT_STATES, c)))
           {
             symptstateB = rbinom(1, *REAL(VECTOR_ELT(PROB_SYMPT, c)));
           }
         }
         if ((symptstateB == 1) & (treatedW == 0))
         {
           treatedW = 1;
         }
         if (treatedW == 1)
         {
           curedW  = rbinom(1, *REAL(VECTOR_ELT(TREATPROBS, location)));
         }
         if (curedW == 1)
         {
           REAL(RES)[i+(*size+11)*j] = 9;  /*recovered: survival health state */
         }
       }
       for (ind=0; ind<j; ind++)
       {
         if (REAL(RES)[i+(*size+11)*(ind)]==location)
         {
           change_state = 0;
         }
       }
       if ((location == 2) & (change_state == 1))
       {
         newCIN1[j] = newCIN1[j] + 1.0;
       }
       if ((location == 3) & (change_state == 1))
       {
         newCIN2[j] = newCIN2[j] + 1.0;
       }
       if ((location == 4) & (change_state == 1))
       {
         newCIN3[j] = newCIN3[j] + 1.0;
       }
       if ((location == 5) & (change_state == 1))
       {
         newCC[j] = newCC[j] + 1.0;
       }
       if (location != REAL(RES)[i+(*size+11)*(j-1)])
       {
         if (location == 1)
         {
           newHPV[j] = newHPV[j] + 1.0;
         }
         if (location == 10)
         {
           newCCMDeath[j] = newCCMDeath[j] + 1.0;
         }
         if (location == 11)
         {
           newOtherDeath[j] = newOtherDeath[j] + 1.0;
         }
       }
       if (treatedW == 1)
       {
        md_cost[i]  = md_cost[i]  + *costCoefs_md;     /* med direct cost of step j for individual i     */
         nmd_cost[i] = nmd_cost[i] + *costCoefs_nmd;    /* non med direct cost of step j for individual i */
         i_cost[i]   = i_cost[i]   + *costCoefs_i;      /* indirect cost of step j for individual i       */
         md_cost_disc[i]  = md_cost_disc[i]  + *costCoefs_md/pow((1+*disc/200), j-1);
         nmd_cost_disc[i] = nmd_cost_disc[i] + *costCoefs_nmd/pow((1+*disc/200), j-1);
         i_cost_disc[i]   = i_cost_disc[i]   + *costCoefs_i/pow((1+*disc/200), j-1);
       }
       /* Set individual indicators */
       if (newCytos[j] > 0)
       {
         CytoW = CytoW+global_scrA.nCyto+global_scrB.nCyto+global_scrC.nCyto+global_scrD.nCyto;
       }
       if (newHPVTests[j] > 0)
       {
         HPVTestW = HPVTestW+global_scrB.nHPV+global_scrC.nHPV+global_scrD.nHPV;
       }
       if (newColpo[j] > 0)
       {
         ColpoW = ColpoW+global_scrA.nColpo+global_scrB.nColpo+global_scrC.nColpo+global_scrD.nColpo;
       }
       if (newBiops[j] > 0)
       {
         BiopW = BiopW+global_scrB.nBiops+global_scrC.nBiops+global_scrD.nBiops;
       }
       if (treatedW != 0)
       {
         TreatedW2 = TreatedW2+1;
       }
       if (symptstateB != 0)
       {
         SymptW = SymptW+1;
       }
     }

     REAL(RES)[j*(*size+10)+j-11] = newHPV[j];                    /* new HPV cases                      */
     REAL(RES)[j*(*size+10)+j-10] = newCC[j];                     /* new CC cases                       */
     REAL(RES)[j*(*size+10)+j-9]  = newCIN1[j];                   /* new CIN1 cases                     */
     REAL(RES)[j*(*size+10)+j-8]  = newCIN2[j];                   /* new CIN2 cases                     */
     REAL(RES)[j*(*size+10)+j-7]  = newCIN3[j];                   /* new CIN3 cases                     */
     REAL(RES)[j*(*size+10)+j-6]  = newCCMDeath[j];               /* new CC death                       */
     REAL(RES)[j*(*size+10)+j-5]  = newOtherDeath[j];             /* new death by other causes          */
     REAL(RES)[j*(*size+10)+j-4]  = newCytos[j];                  /* number of cytologies conducted     */
     REAL(RES)[j*(*size+10)+j-3]  = newColpo[j];                  /* number of colpos conducted         */
     REAL(RES)[j*(*size+10)+j-2]  = newHPVTests[j];               /* number of HPV tests conducted      */
     REAL(RES)[j*(*size+10)+j-1]  = newBiops[j];                  /* number of biopsies conducted       */
    }
    REAL(RES)[i+(*size+11)*(*steps)]    = rbinom(1, *p_men);      /* set sex                            */
    REAL(RES)[i+(*size+11)*(*steps+1)]  = utilities[i];           /* set utilities                      */
    REAL(RES)[i+(*size+11)*(*steps+2)]  = md_cost[i];             /* set direct medical costs           */
    REAL(RES)[i+(*size+11)*(*steps+3)]  = nmd_cost[i];            /* set direct non medical costs       */
    REAL(RES)[i+(*size+11)*(*steps+4)]  = i_cost[i];              /* set indirect costs                 */

    REAL(RES)[i+(*size+11)*(*steps+5)]  = utilities_disc[i];      /* set discounted utilities           */
    REAL(RES)[i+(*size+11)*(*steps+6)]  = md_cost_disc[i];        /* set discounted dir med costs       */
    REAL(RES)[i+(*size+11)*(*steps+7)]  = nmd_cost_disc[i];       /* set discounted dir non med costs   */
    REAL(RES)[i+(*size+11)*(*steps+8)]  = i_cost_disc[i];         /* set discounted indirect costs      */
    REAL(RES)[i+(*size+11)*(*steps+9)]  = le[i];                  /* set life expectancy                */
    REAL(RES)[i+(*size+11)*(*steps+10)] = le_disc[i];             /* set discounted life expectancy     */
    REAL(RES)[i+(*size+11)*(*steps+11)] = vaccinatedW;            /* vaccinated woman                   */
    REAL(RES)[i+(*size+11)*(*steps+12)] = TreatedW2;              /* treated woman                      */
    REAL(RES)[i+(*size+11)*(*steps+13)] = SymptW;                 /* symptoms                           */
    REAL(RES)[i+(*size+11)*(*steps+14)] = CytoW;                  /* Cytology                           */
    REAL(RES)[i+(*size+11)*(*steps+15)] = HPVTestW;               /* HPV test                           */
    REAL(RES)[i+(*size+11)*(*steps+16)] = ColpoW;                 /* Colposcopy                         */
    REAL(RES)[i+(*size+11)*(*steps+17)] = BiopW;                  /* Biopsy                             */
  }

  PutRNGstate();

  UNPROTECT(5);

  return RES;
}
