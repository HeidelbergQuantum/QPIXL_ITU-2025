from helper import *
from qiskit import QuantumCircuit


def cFRQIangs(a, compression, pre_pattern=None,post_pattern=None):
    """
    VARIANT FUNCTION - more of an exampole.
    Constructs the QPIXL circuit.
    This function takes an input image represented as a grayscale array, compresses it
    based on the specified compression ratio, and generates a quantum circuit that 
    encodes the image. The circuit can optionally include pre- and post-pattern operations around controlled rotations.
    Args:
        a (numpy.ndarray): Grayscale image represented as a 1D array of pixel intensities.
        compression (float): Compression ratio as a percentage (0-100). Determines the 
                                proportion of smallest absolute values in the image data 
                                to set to zero.
        pre_pattern (callable, optional): A function that applies a custom operation 
                                            to the circuit before a controlled rotation. 
                                            Defaults to None.
        post_pattern (callable, optional): A function that applies a custom operation 
                                            to the circuit after a controlled rotation. 
                                            Defaults to None.
    Returns:
        QuantumCircuit: A quantum circuit implementing the compressed FRQI representation 
                        of the input image.
    Notes:
        - The input array is first converted to angles and preprocessed.
        - The Walsh-Hadamard Transform (WHT) is applied to the data, followed by 
            a Gray code permutation.
        - The smallest absolute values in the transformed data are set to zero based on 
            the compression parameter.
        """
    a = convertToAngles(a) # convert grayscale to angles
    a = preprocess_image(a) # need to flatten the transpose for easier decoding, 
                            # only really necessary if you want to recover an image.
                            # for classification tasks etc. transpose isn't required.
    n = len(a)
    k = ilog2(n)

    a = 2*a 
    a = sfwht(a)
    a = grayPermutation(a) 
    a_sort_ind = np.argsort(np.abs(a))

    # set smallest absolute values of a to zero according to compression param
    cutoff = int((compression / 100.0) * n)
    for it in a_sort_ind[:cutoff]:
        a[it] = 0
    # print(a)
    # Construct FRQI circuit
    circuit = QuantumCircuit(k + 2)
    # Hadamard register
    circuit.h(range(2,k+2))
    circuit.x(0)
    # Compressed uniformly controlled rotation register
    ctrl, pc, i = 0, 0, 0
    while i < (2**k):
        # Reset the parity check
        pc = int(0)

        # Add RY gate
        if a[i] != 0:
            if pre_pattern is None:
                circuit.ry(a[i],1)
            else:
                pre_pattern(circuit)
                circuit.cry(a[i],0,1)
                post_pattern(circuit)
            
        # Loop over sequence of consecutive zero angles to 
        # cancel out CNOTS (or rather, to not include them)
        if i == ((2**k) - 1):
            ctrl=0
        else:
            ctrl = grayCode(i) ^ grayCode(i+1)
            ctrl = k - countr_zero(ctrl, n_bits=k+1) - 1

        # Update parity check
        pc ^= (2**ctrl)
        i += 1
        
        while i < (2**k) and a[i] == 0:
            # Compute control qubit
            if i == ((2**k) - 1):
                ctrl=0
            else:
                ctrl = grayCode(i) ^ grayCode(i+1)
                ctrl = k - countr_zero(ctrl, n_bits=k+1) - 1

            # Update parity check
            pc ^= (2**ctrl)
            i += 1
                        
        for j in range(k):
            if (pc >> j)  &  1:
                circuit.cx(k-j+1, 1)
    return circuit

def decodeAngQPIXL(state, qc, trace ,max_pixel_val=255, min_pixel_val=0):
    """Automatically decodes qpixl output statevector - taking into account 
    that there are other qubits in the cirucit - this version only returns one image. 

    Args:
        state (statevector array): statevector from simulator - beware of bit ordering
        qc (qiskit circuit): the circuit used for the state generation
        max/min_pixel_val (int, optional): normalization values. Defaults to 255/0.
    Returns:
        np.array: your image, flat
    """
    decoded_data = []
    datum=0
    to_trace = list(range(trace))
    to_trace.pop(trace-datum-1)
    test = decodeQPIXL(partial_trace(state, [qc.qubits.index(qubit) for qubit in [qc.qubits[qub] for qub in to_trace]]).probabilities())
    return convertToGrayscale(np.array([test[permute_bits(i,len(qc.qubits)-trace,datum)] for i in range(len(test))]),max_pixel_val,min_pixel_val)




# def cFRQI(a, compression):
#     """    Takes a standard image in a numpy array (so that the matrix looks like
#     the image you want if you picture the pixels) and returns the QPIXL
#     compressed FRQI circuit. The compression ratio determines
#     how many gates will be filtered and then cancelled out. Made into code from this paper:
#     https://www.nature.com/articles/s41598-022-11024-y

#     Args:
#         a (np.array): numpy array of image, must be flattened and padded with zeros up to a power of two
#         compression (float): number between 0 an 100, where 0 is no compression and 100 is no image

#     Returns:
#         QuantumCircuit: qiskit circuit that prepared the encoded image
#     """
#     a = convertToAngles(a) # convert grayscale to angles
#     a = preprocess_image(a) # need to flatten the transpose for easier decoding, 
#                             # only really necessary if you want to recover an image.
#                             # for classification tasks etc. transpose isn't required.
#     n = len(a)
#     k = ilog2(n)

#     a = 2*a 
#     a = sfwht(a)
#     a = grayPermutation(a) 
#     a_sort_ind = np.argsort(np.abs(a))

#     # set smallest absolute values of a to zero according to compression param
#     cutoff = int((compression / 100.0) * n)
#     for it in a_sort_ind[:cutoff]:
#         a[it] = 0
#     # print(a)
#     # Construct FRQI circuit
#     circuit = QuantumCircuit(k + 1)
#     # Hadamard register
#     circuit.h(range(k))
#     # Compressed uniformly controlled rotation register
#     ctrl, pc, i = 0, 0, 0
#     while i < (2**k):
#         # Reset the parity check
#         pc = int(0)

#         # Add RY gate
#         if a[i] != 0:
#             circuit.ry(a[i], k)

#         # Loop over sequence of consecutive zero angles to 
#         # cancel out CNOTS (or rather, to not include them)
#         if i == ((2**k) - 1):
#             ctrl=0
#         else:
#             ctrl = grayCode(i) ^ grayCode(i+1)
#             ctrl = k - countr_zero(ctrl, n_bits=k+1) - 1

#         # Update parity check
#         pc ^= (2**ctrl)
#         i += 1
        
#         while i < (2**k) and a[i] == 0:
#             # Compute control qubit
#             if i == ((2**k) - 1):
#                 ctrl=0
#             else:
#                 ctrl = grayCode(i) ^ grayCode(i+1)
#                 ctrl = k - countr_zero(ctrl, n_bits=k+1) - 1

#             # Update parity check
#             pc ^= (2**ctrl)
#             i += 1
                        
#         for j in range(k):
#             if (pc >> j)  &  1:
#                 circuit.cx(j, k)
#     return circuit.reverse_bits()
