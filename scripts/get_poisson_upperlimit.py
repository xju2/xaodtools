#!/usr/bin/env python
"""
get upper limit from simple Poisson distribution
use high mass workspace as input
"""

import ROOT
import sys
def get_limit(num_expected, cls=0.05):
    limit = num_expected
    if num_expected == 0:
        step = 0.1
    else:
        step = num_expected * 0.1

    cdf = ROOT.Math.poisson_cdf_c(limit, num_expected)
    n_tries = 0
    while cdf > cls and n_tries < 100:
        limit += step
        cdf = ROOT.Math.poisson_cdf_c(limit, num_expected)
        n_tries += 1
        #print cdf,limit

    print "CL",1-cls,"% upper limit:",limit," for",num_expected
    return limit

def process(file_name, mH, lumi=13.3):
    ws_name = "combined"
    mc_name = "ModelConfig"
    pdf_name = "simPdf"
    data_name = "obsData"
    threshold =  0.68

    prod_name = "ggF"
    pdf_4mu_name = "ATLAS_Signal_"+prod_name+"_ggF_4mu_13TeV_cbga"

    fin = ROOT.TFile.Open(file_name)
    ws = fin.Get(ws_name)
    pdf = ws.obj(pdf_name)
    observables = ws.obj(mc_name).GetObservables()
    obs = observables.first()
    data = ws.obj(data_name)

    ws.var("mH").setVal(mH)
    out = ""
    out += "mH: "+str(mH)+"\n"

    ## obtain acceptance
    acceptance = 0.
    categories = pdf.indexCat()
    cat_iter = ROOT.TIter(categories.typeIterator())
    obj = cat_iter()
    while obj:
        cat_name = obj.GetName()
        poly_name = "ACpol_"+prod_name+"_"+cat_name.replace("Cat","")
        var = ws.obj(poly_name)
        if not var:
            print poly_name,"does not exist"
            obj = cat_iter()
            continue

        acceptance += var.getVal()
        obj = cat_iter()

    out += "acceptance: "+str(acceptance)+"\n"

    ## find 95% mass window for mH
    width = mH*0.02
    step = width*0.01
    ## take 4mu channel as reference, as it has worst resolution
    pdf_4mu = ws.obj(pdf_4mu_name)
    pdf_4mu_total = pdf_4mu.createIntegral(observables).getVal()
    #print "4mu total:",pdf_4mu_total
    obs.setRange("res0", mH-width, mH+width)
    frac = pdf_4mu.createIntegral(observables, ROOT.RooFit.Range("res0")).getVal()/pdf_4mu_total
    n_tries = 0
    low_mass = -1
    high_mass = 1E9
    range_def = ""

    while frac < threshold and n_tries < 100:
        n_tries += 1
        range_name = "res"+str(n_tries)
        range_def = range_name
        low_mass = mH - width - step*n_tries
        high_mass = mH + width + step*n_tries
        obs.setRange(range_name, low_mass, high_mass)
        frac = pdf_4mu.createIntegral(observables, ROOT.RooFit.Range(range_name)).getVal()/pdf_4mu_total

    out += "mass range: {:.1f},{:.1f} with fraction: {:.3f} and threshold: {:.2f}\n".format(low_mass,high_mass,frac,threshold)

    ## find number of observed events
    datasets = data.split(categories, True)
    num_obs = 0
    for data_ch in datasets:
        for ientry in range(data_ch.numEntries()):
            m4l = data_ch.get(ientry).first().getVal()
            if m4l < low_mass or m4l > high_mass:
                continue
            num_obs += 1

    out += "observed events:"+str(num_obs)+"\n"

    # find expected number of background
    obs = ws.obj(mc_name).GetParametersOfInterest()
    obs_itr = ROOT.TIter(obs.createIterator())
    obj = obs_itr()
    while obj:
        obj.setVal(0.)
        obj = obs_itr()

    cat_iter.Reset()
    obj = cat_iter()
    num_total_bkg = 0
    while obj:
        pdf_ch = pdf.getPdf(obj.GetName())
        num_total_bkg += pdf_ch.expectedEvents(observables)
        obj = cat_iter()


    total_bkg_integral = pdf.createIntegral(observables).getVal()
    #num_total_bkg = pdf.expectedEvents(observables)
    print "total background:", num_total_bkg
    frac_bkg = pdf.createIntegral(observables, ROOT.RooFit.Range(range_def)).getVal()
    num_bkg = frac_bkg/total_bkg_integral*num_total_bkg
    out += "number of background: {:.2f}\n".format(num_bkg)

    # get upper limit
    poisson_limit = get_limit(num_obs)
    xs_limit = (poisson_limit-num_bkg)/acceptance/lumi
    out += "XS limit: "+str(xs_limit)+"\n"

    print out

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print sys.argv[0],"arg "
        print "arg: get_limit number_expected cls;"
        print "     process file_name mH lumi"
        sys.exit(1)

    arg = sys.argv[1]
    if arg == "get_limit":
        if len(sys.argv) < 3:
            print "get_limit number_expected cls"
            sys.exit(2)

        num_exp = int(sys.argv[2])
        if len(sys.argv) > 3:
            cls = float(sys.argv[3])
        else:
            cls = 0.05

        get_limit(num_exp, cls)
    elif arg == "process":
        if len(sys.argv) < 4:
            print "process file_name mH (lumi)"
            sys.exit(3)

        file_name = sys.argv[2]
        mH = float(sys.argv[3])
        if len(sys.argv) > 4:
            lumi = float(sys.argv[4])
        else:
            lumi = 13.3

        process(file_name, mH, lumi)
    else:
        print "I don't know"

