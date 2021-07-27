"""
this main kicks off the preprocessing for the vizgen data.
It prepares and shapes the data in the appropriate form so they can be consumed by the celltyping algo
The output of this function is the input to pciSeq
"""
import logging
from pciSeq.src.preprocess.vizgen import get_data
logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

if __name__ == "__main__":
    merfish_id = 'MERFISH_M_Z'
    slice_ids = [
        "MsBrain_ZM0_VS6_JH_V6_05-15-2021",
        "MsBrain_ZM1_VS6_JH_V11_05-16-2021",
        "MsBrain_ZM2_VS6_JH_V11_05-15-2021",
        "MsBrain_ZM3_VS6_JH_V11_05-17-2021",
        "MsBrain_ZM4_VS6_JH_V11_05-11-2021",
        "MsBrain_ZM5.1_VS6_JH_V11_05-12-2021",
        "MsBrain_ZM5.2_VS6_JH_V6_05-13-2021",
        "MsBrain_ZM6.1_VS6_V6_JH_05-11-2021",
        "MsBrain_ZM7.1_VS6_V6_JH_05-12-2021",
        "MsBrain_ZM7.2_VS6_JH_V11_05-13-2021"
        ]
    region_ids = ['region_0']

    for slice_id in slice_ids:
        for region_id in region_ids:
            logger.info("\n Started slice %s, region %s" % (slice_id, region_id))
            get_data.run(merfish_id, slice_id, region_id)


    logger.info('Done!')

