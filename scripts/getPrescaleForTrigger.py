#!/usr/bin/env python
#asetup.sh 20.2.3,here
import sys
from collections import defaultdict
from PyCool import cool
import TriggerMenu.l1.PrescaleHelper as psh

import plot_data_per_run

dbSvc = cool.DatabaseSvcFactory.databaseService()

debug =False
def getCoolData(run, db, folder, fields):

    db = dbSvc.openDatabase('oracle://ATLAS_COOLPROD;schema=%s;dbname=CONDBR2' % db,False)
    fld = db.getFolder( folder )
    objs = fld.browseObjects( (run<<32)+1, ((run+1)<<32)-1, cool.ChannelSelection() )

    result = []
    while objs.goToNext():
        obj = objs.currentRef()
        result += [ (obj.since() & 0xffffffff, (obj.until()-1) & 0xffffffff, obj.channelId(), map(obj.payloadValue, fields)) ]

    return result


def getChainsInfo(run, triggers):

    hltmenu = getCoolData( run,
                           'ATLAS_COOLONL_TRIGGER',
                           '/TRIGGER/HLT/Menu',
                           ["ChainName","LowerChainName","ChainCounter"] )

    return dict([ (v[0],v) for (b,e,c,v) in hltmenu if v[0] in triggers])


def getItemsInfo(run, triggers):

    lvl1menu = getCoolData( run,
                            'ATLAS_COOLONL_TRIGGER',
                            '/TRIGGER/LVL1/Menu',
                            ["ItemName"] )

    return dict([(v[0],[c]+v) for (b,e,c,v) in lvl1menu if v[0] in triggers])


def getL1PSK(run):

    return getCoolData( run,
                        'ATLAS_COOLONL_TRIGGER',
                        '/TRIGGER/LVL1/Lvl1ConfigKey',
                        ['Lvl1PrescaleConfigurationKey'] )


def getHLTPSK(run):

    return getCoolData( run,
                        'ATLAS_COOLONL_TRIGGER',
                        '/TRIGGER/HLT/PrescaleKey',
                        ['HltPrescaleKey'] )


def getL1Prescales(run, ctpIds):
    ps = getCoolData( run,
                      'ATLAS_COOLONL_TRIGGER',
                      '/TRIGGER/LVL1/Prescales',
                      ['Lvl1Prescale'] )

    psd = defaultdict(list)
    [map(psd[c].append,((b,e,int(v[0])),) ) for (b,e,c,v) in ps if c in ctpIds]
    return psd
    

def getHLTPrescales(run, chainIds):
    ps = getCoolData( run,
                      'ATLAS_COOLONL_TRIGGER',
                      '/TRIGGER/HLT/Prescales',
                      ['Prescale'] )

    psd = defaultdict(list)
    [map(psd[c-20000].append,((b,e,float(v[0])),) ) for (b,e,c,v) in ps if (c-20000) in chainIds]
    return psd


def getCombinedPrescales(triggerIDMap, l1pss, hltpss):

    psmap = {}

    for trigger, idmap in sorted(triggerIDMap.items()):
        if debug: print '\n',trigger
        if debug: print '-' * len(trigger)

        l1ID  = idmap["l1Id"]
        hltID = idmap["hltId"]

        l1lbps = l1pss[l1ID]
        hltlbps = hltpss[hltID]

        lbstart = sorted(list(set([v[0] for v in l1lbps]).union([v[0] for v in hltlbps])))
        lbstart.remove(0)

        psmap[trigger] = []
        for lb in lbstart:

            l1ps = None
            hltps = None

            for b,e,v in l1lbps:
                if lb>=b and lb<=e:
                    l1ps = psh.getPrescaleFromCut(v)
                    break
            for b,e,v in hltlbps:
                if lb>=b and lb<=e:
                    hltps = v
                    break

            ps = abs(l1ps * hltps)
            if l1ps<0 or hltps<0:
                ps = -1 * ps

            psmap[trigger] += [ (lb, ps, l1ps, hltps) ]

            if debug: print "LB %3i (l1ps %10.1f  hltps %7.2f)  ==>  ps %10.1f" % (lb,l1ps,hltps,ps)

    return psmap




def main(run=284484):

    jet_triggers = [
        "HLT_j15",
        "HLT_j25",
        "HLT_j60",
        "HLT_j100",
        "HLT_j110",
        "HLT_j150",
        "HLT_j175",
        "HLT_j200",
        "HLT_j260",
        "HLT_j300",
        "HLT_j320",
        "HLT_j360",
        "HLT_j400",
        "HLT_j420",
    ]
    muon_triggers = [
        # Di muon
        "HLT_2mu10", "HLT_mu18_mu8noL1", "HLT_2mu10_nomucomb",
        "HLT_mu20_mu8noL1", "HLT_mu20_nomucomb_mu6noL1_nscan03",
        "HLT_2mu14_nomucomb", "HLT_mu22_mu8noL1", "HLT_2mu14",
        # Tri-muons
        "HLT_3mu6", "HLT_3mu6_msonly", "HLT_mu18_2mu4noL1",
        "HLT_mu20_2mu4noL1", "HLT_3mu4", "HLT_mu6_2mu4",
        "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6",
        "HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb",
        "HLT_3mu6", "HLT_3mu6_msonly"
    ]

    #triggers = jet_triggers
    triggers = muon_triggers

    chains = getChainsInfo( run, triggers )

    noL1Input = [c for c,v in chains.items() if v[1]=='']
    if noL1Input:
        print "These chains have no L1 input (ie they run on all L1Accept) and should not be used. They are removed from further analysis: %s" % ", ".join(noL1Input)
        map(triggers.remove, noL1Input)

    items = getItemsInfo( run, [v[1] for v in chains.values()] )

    triggerIDMap = dict([ (c, {"l1Id" : int(items[chains[c][1]][0]) , "hltId" : int(chains[c][2])}) for c in chains])

    # item and chain IDs (for getting the correct prescales)
    itemIDs  = [idd["l1Id"] for idd in triggerIDMap.values()]
    chainIDs = [idd["hltId"] for idd in triggerIDMap.values()]

    # getting the prescales
    l1pss  = getL1Prescales( run, itemIDs )
    hltpss = getHLTPrescales( run, chainIDs )

    psmaps = getCombinedPrescales(triggerIDMap, l1pss, hltpss)
    out = ""
    for key,value in psmaps.iteritems():
        out +='['+ key+']\n'
        for each_lb in value:
            lb,ps,l1ps,hltps = each_lb
            out += "LB %3i, ps %10.1f" % (lb,ps) + "\n"

    with open(str(run)+"_ps.txt", 'w') as f:
        f.write(out)


if __name__=="__main__":
    if len(sys.argv) < 2:
        lumi_dic = plot_data_per_run.load_lumi()
        print lumi_dic
        for key,value in lumi_dic.iteritems():
            main(int(key))
    else:
        main(int(sys.argv[1]))
