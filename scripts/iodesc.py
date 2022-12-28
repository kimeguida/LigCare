

# Copyright (c) 2022 LIT
# For academic use only.
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import numpy as np



def write_desc(indexes_, fingerprints_, ofile_):
    descriptors = []
    for i in indexes_:
        concat_fp = [i]
        for fp in fingerprints_:
            concat_fp += list(fp[i])
        concat_fp = np.array(concat_fp)
        descriptors.append(concat_fp)
    np.save(ofile_, descriptors)




def read_desc(desc_, pos_=None):
    data = np.load(desc_)
    indexes = data[:, :1]
    indexes = [int(i[0]) for i in indexes]
    descriptors = data[:, 1:]

    if pos_ is not None:
        extracted_columns = []
        for p_inf, p_sup in pos_:
            extracted_columns.append(descriptors[:, list(range(p_inf, p_sup+1))])
            # implemented as such to handle [x1, x1] col extraction whithout numpy repeating 
        #print(extracted_columns)

        if len(pos_) == 1:
            reduced_descriptors = extracted_columns[0]
        else:
            reduced_descriptors = np.concatenate(extracted_columns, axis=1)

        descriptors = reduced_descriptors


    return indexes, descriptors




def write_bsa(bsa_, ofile_):
    with open(ofile_, 'w') as of:
        of.write('mol2_index\tbsa\n')
        for i, bsa in bsa_.items():
            of.write(f'{i+1}\t{bsa}\n')





def read_bsa(file_):
    buriedness = {}
    with open(file_, 'r') as f:
        for l in f:
            if l.startswith('mol2_index'):
                continue
            cols = l.split('\t')
            idx = int(cols[0]) - 1
            bsa = int(cols[1])
            buriedness[idx] = bsa
    return buriedness







def write_labels(cav_labels_, ofile_):
    np.save(ofile_, cav_labels_)




def read_labels(cav_labels_):
    cav_labels = np.load(cav_labels_)
    return cav_labels

