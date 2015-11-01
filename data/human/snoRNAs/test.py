import sys
sys.path.append("../../../modules")
import snoRNA

cd_snors = snoRNA.read_snoRNAs_from_table(sys.argv[1], None, True)
