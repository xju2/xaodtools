#!/usr/bin/env python
"""
Get yields and shape for the reducible background in H4l analysis
"""

import ROOT
from optparse import OptionParser


def run(options):
    ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
    from ROOT import H4lBackgroundReader

    b = H4lBackgroundReader()
    b.loadShapes("$ROOTCOREBIN/../H4lBackgroundReader/data/shapes_SignalRegion.root");
    b.loadNorms("$ROOTCOREBIN/../H4lBackgroundReader/data/normalizations_m4lFullRange_1ifb.root");
    b.setLumiInvFb(36.1);

    if options.no_VBF:
        channels = [ROOT.H4lBkg.Channel._4l, ROOT.H4lBkg.Channel._4e, ROOT.H4lBkg.Channel._2mu2e, ROOT.H4lBkg.Channel._4mu, ROOT.H4lBkg.Channel._2e2mu]
        for channel in channels:
            bgm = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, ROOT.H4lBkg.Category.prodInclusive,  channel, 130, 1500);
            m = bgm.getNorm()
            print ROOT.H4lBkg.Channel.toString( channel ), '\t', m.getString()

    else:
        categories = [ ROOT.H4lBkg.Category._HMggF, ROOT.H4lBkg.Category._HMVBF ]
        for category in categories:
            print '-----------',ROOT.H4lBkg.Category.toString(category),'------------'
            bgm_4mu = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, category,  ROOT.H4lBkg.Channel._4mu , 130, 1500);
            bgm_4e = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, category,  ROOT.H4lBkg.Channel._4e , 130, 1500);
            bgm_2e2mu = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, category,  ROOT.H4lBkg.Channel._2e2mu , 130, 1500);
            bgm_2mu2e = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, category,  ROOT.H4lBkg.Channel._2mu2e , 130, 1500);
            bgm_4l = b.getBkgMeasurement( ROOT.H4lBkg.Observable.m4lHM, category,  ROOT.H4lBkg.Channel._4l, 130, 1500);
            bgm_mixed = ROOT.Add(bgm_2mu2e,bgm_2e2mu, category, ROOT.H4lBkg.Channel._4l);

            print ROOT.H4lBkg.Channel.toString( ROOT.H4lBkg.Channel._4mu ), '\t', bgm_4mu.getNorm().getString()
            print ROOT.H4lBkg.Channel.toString( ROOT.H4lBkg.Channel._4e ), '\t', bgm_4e.getNorm().getString()
            print '2e2mu', '\t', bgm_mixed.getNorm().getString()
            print ROOT.H4lBkg.Channel.toString( ROOT.H4lBkg.Channel._4l ), '\t', bgm_4l.getNorm().getString()

            ## example of getting shape
            ## h_m4mu_shape = bgm_4mu.getHistShape();


if __name__ == "__main__":
    usage = "%prog [options]"
    version="%prog 1.1"
    parser = OptionParser(usage=usage, description="get yields for WS", version=version)
    parser.add_option("--noVBF", dest="no_VBF", default=False, action='store_true', help='For Highmass, no VBF-like category')
    (options,args) = parser.parse_args()

    run(options)
