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
    merfish_id = 'MERFISH_F_F'
    slice_ids = [
        "VS6MsBrain_F5_VS6libary_V10_LH_03-17-21"
        ]
    region_ids = ['region_0', 'region_1']

    for slice_id in slice_ids:
        for region_id in region_ids:
            logger.info("\n Started slice %s, region %s" % (slice_id, region_id))
            get_data.run(merfish_id, slice_id, region_id)


    logger.info('Done!')

