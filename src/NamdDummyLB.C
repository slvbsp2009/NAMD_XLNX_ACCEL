
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
#include "NamdDummyLB.h"
#include "NamdDummyLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

void CreateNamdDummyLB() {
  int seqno = LdbInfra::Object()->getLoadbalancerTicket();
  loadbalancer = CProxy_NamdDummyLB::ckNew(CkLBOptions(seqno));
}

NamdDummyLB *AllocateNamdDummyLB() {
  return new NamdDummyLB((CkMigrateMessage*)NULL);
}

NamdDummyLB::NamdDummyLB(CkMigrateMessage *msg): CentralLB(msg) {
  lbname = (char*)"NamdDummyLB";
}
 
NamdDummyLB::NamdDummyLB(const CkLBOptions& opt): CentralLB(opt) {
  lbname = (char*)"NamdDummyLB";
  if (CkMyPe() == 0)
    CkPrintf("[%d] DummyLB created\n",CkMyPe());
}

bool NamdDummyLB::QueryBalanceNow(int _step) {
  return true;
}

bool NamdDummyLB::QueryDumpData() {
  return false;
}

// Dummy work function

void NamdDummyLB::work(LDStats* stats) {
  // CkPrintf("[%d] NamdDummyLB At WORK\n",CkMyPe());
}
