from ligandparam.parametrization import Recipe
from ligandparam.stages import *
#
# class RNALigand(Recipe):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         return
#     def setup(self):
#         raise NotImplementedError("This class is not yet implemented.")
#         self.stages = [
#             StageInitialize("Initialize",,,
#             """
#             StageNormalizeCharge("Normalize1", base_cls=self,
#                     in_mol2=self.name+".antechamber.mol2",
#                     out_mol2=self.name+".antechamber.mol2"),
#             StageGaussian("Minimize", base_cls=self),
#             StageLazyResp("LazyResp", base_cls=self),
#             StageNormalizeCharge("Normalize2", base_cls=self,
#                     in_mol2=self.name+".resp.mol2",
#                     out_mol2=self.name+".resp.mol2"),
#             StageUpdate("UpdateNames", base_cls=self,
#                     in_mol2=self.name+'.antechamber.mol2',
#                     to_update=self.name+'.resp.mol2',
#                     out_mol2=self.name+'.resp.mol2',
#                     update_names=True,
#                     update_types=False,
#                     update_resname=True),
#             StageParmChk("ParmChk", base_cls=self),
#             StageLeap("Leap", base_cls=self)
#             """
#         ]