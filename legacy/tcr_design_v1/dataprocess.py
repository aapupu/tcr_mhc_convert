import os
import pickle
import re

import numpy as np
import pandas as pd

def mk_aa_dict():
    """Returns the amino acid dictionary mapping."""
    amino_acid_dict = {
        '[PAD]': 0,
        'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5,
        'Q': 6, 'E': 7, 'G': 8, 'H': 9, 'I': 10,
        'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
        'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20,
        'X':21, '[MASK]': 22, '[CLS]': 23, '[SEP]': 24,
        '[UNK]':25
    }
    return amino_acid_dict

def mk_bv_dict():
    bv_dict = {
        '[PAD]': 0,'TRBV1': 1,'TRBV10-1': 2,'TRBV10-2': 3,'TRBV10-3': 4,'TRBV11-1': 5,'TRBV11-2': 6,'TRBV11-3': 7,'TRBV12-1': 8,'TRBV12-2': 9,
        'TRBV12-3': 10,'TRBV12-4': 11,'TRBV12-5': 12,'TRBV13': 13,'TRBV14': 14,'TRBV15': 15,'TRBV16': 16,'TRBV17': 17,'TRBV18': 18,'TRBV19': 19,
        'TRBV2': 20,'TRBV20-1': 21,'TRBV20/OR9-2': 22,'TRBV21-1': 23,'TRBV21/OR9-2': 24,'TRBV22-1': 25,'TRBV22/OR9-2': 26,'TRBV23-1': 27,'TRBV23/OR9-2': 28,
        'TRBV24-1': 29,'TRBV24/OR9-2': 30,'TRBV25-1': 31,'TRBV26': 32,'TRBV26/OR9-2': 33,'TRBV27': 34,'TRBV28': 35,'TRBV29-1': 36,'TRBV29/OR9-2': 37,'TRBV3-1': 38,
        'TRBV30': 39,'TRBV4-1': 40,'TRBV4-2': 41,'TRBV4-3': 42,'TRBV5-1': 43,'TRBV5-2': 44,'TRBV5-3': 45,'TRBV5-4': 46,'TRBV5-5': 47,'TRBV5-6': 48,'TRBV5-7': 49,
        'TRBV5-8': 50,'TRBV6-1': 51,'TRBV6-2': 52,'TRBV6-4': 53,'TRBV6-5': 54,'TRBV6-6': 55,'TRBV6-7': 56,'TRBV6-8': 57,'TRBV6-9': 58,'TRBV7-1': 59,'TRBV7-2': 60,
        'TRBV7-3': 61,'TRBV7-4': 62,'TRBV7-5': 63,'TRBV7-6': 64,'TRBV7-7': 65,'TRBV7-8': 66,'TRBV7-9': 67,'TRBV8-1': 68,'TRBV8-2': 69,'TRBV9': 70,'TRBVA': 71,
        'TRBVB': 72, 'TRBV6-3': 73, 'TRBV3-2':74, 'TRBVA/OR9-2':75, 'TRBV25/OR9-2':76, '[UNK]':77
    }
    return bv_dict

def mk_bj_dict():
    aj_dict = {
        '[PAD]': 0, 'TRBJ1-1': 1,'TRBJ1-2': 2,'TRBJ1-3': 3,'TRBJ1-4': 4,'TRBJ1-5': 5,'TRBJ1-6': 6,
        'TRBJ2-1': 7,'TRBJ2-2': 8,'TRBJ2-2P': 9,'TRBJ2-3': 10,
        'TRBJ2-4': 11,'TRBJ2-5': 12,'TRBJ2-6': 13,'TRBJ2-7': 14, '[UNK]':15
    }
    return aj_dict

def mk_av_dict():
    av_dict = {
        '[PAD]': 0,'TRAV1-1': 1,'TRAV1-2': 2,'TRAV10': 3,'TRAV11': 4,'TRAV12-1': 5,'TRAV12-2': 6,'TRAV12-3': 7,'TRAV13-1': 8,'TRAV13-2': 9,
        'TRAV14/DV4': 10,'TRAV15': 11,'TRAV16': 12,'TRAV17': 13,'TRAV18': 14,'TRAV19': 15,'TRAV2': 16,'TRAV20': 17,'TRAV21': 18,
        'TRAV22': 19,'TRAV23/DV6': 20,'TRAV24': 21,'TRAV25': 22,'TRAV26-1': 23,'TRAV26-2': 24,'TRAV27': 25,'TRAV28': 26,
        'TRAV29/DV5': 27,'TRAV3': 28,'TRAV30': 29,'TRAV31': 30,'TRAV32': 31,'TRAV33': 32,'TRAV34': 33,'TRAV35': 34,
        'TRAV36/DV7': 35,'TRAV37': 36,'TRAV38-1': 37,'TRAV38-2/DV8': 38,'TRAV39': 39,'TRAV4': 40,'TRAV40': 41,
        'TRAV41': 42,'TRAV5': 43,'TRAV6': 44,'TRAV7': 45,'TRAV8-1': 46,'TRAV8-2': 47,'TRAV8-3': 48,
        'TRAV8-4': 49,'TRAV8-5': 50,'TRAV8-6': 51,'TRAV8-7': 52,'TRAV9-1': 53,'TRAV9-2': 54, '[UNK]':55
    }
    return av_dict

def mk_aj_dict():
    aj_dict = {
        '[PAD]': 0,'TRAJ1': 1,'TRAJ10': 2,'TRAJ11': 3,'TRAJ12': 4,'TRAJ13': 5,'TRAJ14': 6,'TRAJ15': 7,'TRAJ16': 8,
        'TRAJ17': 9,'TRAJ18': 10,'TRAJ19': 11,'TRAJ2': 12,'TRAJ20': 13,'TRAJ21': 14,'TRAJ22': 15,'TRAJ23': 16,
        'TRAJ24': 17,'TRAJ25': 18,'TRAJ26': 19,'TRAJ27': 20,'TRAJ28': 21,'TRAJ29': 22,'TRAJ3': 23,'TRAJ30': 24,
        'TRAJ31': 25,'TRAJ32': 26,'TRAJ33': 27,'TRAJ34': 28,'TRAJ35': 29,'TRAJ36': 30,'TRAJ37': 31,'TRAJ38': 32,
        'TRAJ39': 33,'TRAJ4': 34,'TRAJ40': 35,'TRAJ41': 36,'TRAJ42': 37,'TRAJ43': 38,'TRAJ44': 39,'TRAJ45': 40,
        'TRAJ46': 41,'TRAJ47': 42,'TRAJ48': 43,'TRAJ49': 44,'TRAJ5': 45,'TRAJ50': 46,'TRAJ51': 47,'TRAJ52': 48,
        'TRAJ53': 49,'TRAJ54': 50,'TRAJ55': 51,'TRAJ56': 52,'TRAJ57': 53,'TRAJ58': 54,'TRAJ59': 55,'TRAJ6': 56,
        'TRAJ60': 57,'TRAJ61': 58,'TRAJ7': 59,'TRAJ8': 60,'TRAJ9': 61, '[UNK]':62
    }
    return aj_dict

def load_mhc_dict(mhc_type, pseudo=True):
    current_dir = os.path.dirname(__file__)
    project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
    
    assert mhc_type in ('HLAI', 'HLAII'),'mhc_type must be HLAI or HLAII'
    if pseudo:
        file_name = 'pseudo_HLAI.csv' if mhc_type == 'HLAI' else 'pseudo_HLAII.csv'
    else:
        file_name = 'esm_emb_HLAI.pkl' if mhc_type == 'HLAI' else 'esm_emb_HLAII.pkl'
    file_path = os.path.join(project_root, 'doc', file_name)
    if pseudo:
        mhc_data = pd.read_csv(file_path, index_col=0)
        mhc_dict = dict(zip(mhc_data['MHC'], mhc_data['pseudo_seq']))
    else:
        with open(file_path, 'rb') as f:
            mhc_dict = pickle.load(f)
    return mhc_dict

def aa_to_vec(aa_seq, aa_dict):
    """Converts an amino acid sequence to a vector representation."""
    aa_seq = aa_seq.replace(u'\xa0', u'').upper()
    return np.array([aa_dict.get(aa, aa_dict['X']) for aa in aa_seq], dtype=int)

def pad_1d(arr, max_len, pad_value=0, dtype=int):
    arr = np.asarray(arr, dtype=dtype)
    if arr.shape[0] >= max_len:
        return arr[:max_len]
    if arr.shape[0] < max_len:
        return np.pad(arr, (0, max_len - arr.shape[0]),
                      mode='constant', constant_values=pad_value)
    return arr

def pad_2d(mat, target_h, target_w, pad_value=0, dtype=float):
    mat = np.asarray(mat, dtype=dtype)
    if mat.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {mat.shape}")

    h, w = mat.shape

    h_trunc = min(h, target_h)
    w_trunc = min(w, target_w)
    mat_cropped = mat[:h_trunc, :w_trunc]

    pad_h = max(target_h - h_trunc, 0)
    pad_w = max(target_w - w_trunc, 0)

    # --- pad ---
    if pad_h > 0 or pad_w > 0:
        mat_padded = np.pad(
            mat_cropped,
            pad_width=((0, pad_h), (0, pad_w)),
            mode="constant",
            constant_values=pad_value,
        )
    else:
        mat_padded = mat_cropped

    if mat_padded.shape != (target_h, target_w):
        raise RuntimeError(
            f"pad_2d failed: got {mat_padded.shape}, expected ({target_h}, {target_w})"
        )
    return mat_padded
    
def mhc_to_aa(mhc_name, mhc_dict):
    mhc_aa = mhc_dict.get(mhc_name)
    if mhc_aa is None:
        raise KeyError(f"MHC name '{mhc_name}' not found in mhc_dict.")
    return mhc_aa

def mhc_to_esm(mhc_name, esm_dict):
    mhc_esm_emb = esm_dict.get(mhc_name)
    if mhc_esm_emb is None:
        raise KeyError(f"MHC name '{mhc_name}' not found in esm_dict.")
    
    pad_rows = 369 - mhc_esm_emb.shape[0]
    if pad_rows > 0:
        padding = np.zeros((pad_rows, mhc_esm_emb.shape[1]), dtype=mhc_esm_emb.dtype)
        mhc_esm_emb = np.concatenate([mhc_esm_emb, padding], axis=0) 
        
    return mhc_esm_emb

def labelMap(df, label):
    return df[label].to_numpy()


def load_gene_aa(csv_path):
    df = pd.read_csv(csv_path)
    return dict(zip(df["gene_name"], df["amino_acid"]))

def clean_aa(seq):
    import re
    if seq is None:
        return ""
    seq = str(seq).strip().upper()
    seq = seq.replace(".", "").replace("-", "").replace(" ", "")
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYX]", "", seq)
    return seq


def find_overlap(a, b, min_overlap=1, max_overlap=30):
    """
    找 a 的 suffix 和 b 的 prefix 的最大重叠长度
    """
    a = clean_aa(a)
    b = clean_aa(b)

    if not a or not b:
        return 0, ""

    max_k = min(len(a), len(b), max_overlap)

    for k in range(max_k, min_overlap - 1, -1):
        if a[-k:] == b[:k]:
            return k, a[-k:]

    return 0, ""


def choose_trbc_gene(trbj, default="TRBC1"):
    """
    Choose beta constant gene by TRBJ family.

    Human TRBJ1-* segments usually pair with TRBC1, while TRBJ2-* segments
    usually pair with TRBC2. Unknown/irregular names fall back to TRBC1.
    """
    trbj = str(trbj).strip().upper()
    if trbj.startswith("TRBJ2"):
        return "TRBC2"
    if trbj.startswith("TRBJ1"):
        return "TRBC1"
    return default
