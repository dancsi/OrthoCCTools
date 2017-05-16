import sys
from pathlib import Path
from subprocess import run, DEVNULL

from bruteforce_cutoff import bruteforce
from cartesian_product import product_fasta
from reconstruct_orthoset import reconstruct


def check_path_exists(path):
    if not path.exists():
        print("Path", str(path), "does not exist")
        sys.exit(1)


def compute_orthoset(fasta_path, output_name):
    print('Computing orthoset for', fasta_path)

    output_pairs = temp_dir / (output_name + '.pairs')
    score_basepath = output_pairs.with_suffix('')

    run([r'..\..\build-win\fastscore', str(fasta_path),
         '--basename={}'.format(score_basepath)])
    bin_path = score_basepath.with_suffix('.bin')
    print('Calculated', bin_path)

    (binding_cutoff, nonbinding_cutoff) = bruteforce(bin_path, fasta_path)

    run([r'..\..\build-win\solver', str(bin_path), '--binding-cutoff={}'.format(
        binding_cutoff), '--nonbinding-cutoff={}'.format(nonbinding_cutoff), '--out-name={}'.format(output_pairs), '--fasta-name={}'.format(fasta_path)], stdout=DEVNULL, stderr=DEVNULL)

    orthoset_fasta_path = reconstruct(output_pairs, fasta_path=fasta_path)

    return orthoset_fasta_path


if len(sys.argv) < 4:
    print('Iterative set building script')
    print('USAGE: iterative_set_builder.py BASE EXTENSION NUM')
    sys.exit(1)

base_fasta_path = Path(sys.argv[1])
extension_fasta_path = Path(sys.argv[2])
num = int(sys.argv[3])

check_path_exists(base_fasta_path)
check_path_exists(extension_fasta_path)

temp_dir = Path('.') / "temp"
temp_dir.mkdir(exist_ok=True)

current_fasta_path = compute_orthoset(base_fasta_path, "initial")

for i in range(num):
    product_fasta_path = temp_dir / '{}.fasta'.format(i + 1)
    product_fasta([(current_fasta_path, 1),
                   (extension_fasta_path, 1)], product_fasta_path)
    current_fasta_path = compute_orthoset(product_fasta_path, str(i + 1))
