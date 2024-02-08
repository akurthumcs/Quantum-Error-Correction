#= This program simulates the process of quantum error correction for a single qubit, utlizing the
technique of repetition codes. First, we create a random qubit, represented as a vector. We then
entangle the qubit with two other qubits, initialized to state |0>, to create a three qubit system,
where |000> and |111> represent logical |0> and logical |1> respectively. We simulate possible error
by sending the newly encoded qubit through a bit-flip channel, which applies a bit flip to at most one
of the physicial qubits. Given this qubit, after it is sent through the error channel, and three 
projective operators, we are able to calculate the syndrome of the new qubit system, we use a parity
matrix to compute the final syndrome, then given the syndrome, we can find the error, and correct it,
by applying the corresponding bit flip operator, the final result is our original logical qubit. =#

#= Here we have a list of macros, which represent operations that we use or compose operations we use
=#
#Gate in first time step
G_1 = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0]
#Gate in second time step
G_2=[1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 1 0 0; 
0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 1 0]

#Projective Operators for Syndrome Calculation
P_1 = [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0;
0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]

P_2 = [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0;
0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]

P_3 = [0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0]

#This is the parity matrix, represents taking parity measurement of second and third qubits and
#first and third qubits
parity_matrix = [0 1 1; 1 0 1]

#Other Operators Used
Pauli_X = [0 1; 1 0]
I = [1 0; 0 1]

#BitFlip Operations
BitFlip1 = kron(kron(Pauli_X, I), I)
BitFlip2 = kron(kron(I, Pauli_X), I)
BitFlip3 = kron(kron(I, I), Pauli_X)

#Here we have some functions used in the code

#Take the tensor product of two vectors
function tensor(v1 :: Array, v2:: Array):: Array
    w = []
    for i = 1:length(v1)
        for j = 1:length(v2)
            push!(w, (v1[i] * v2[j]))
        end
    end
    return w
end

#Send a qubit through a bit flip error channel with probability p, we make the assumption that each
#qubit is acted on individually
function bitFlipChannel(v:: Array, p):: Array
    if rand() < p
        return BitFlip1 * v
    end
    if rand() < p
        return BitFlip2 * v
    end
    if rand() < p
        return BitFlip3 * v
    end
    return v
end
#=Calculate syndrome using projection operators, we apply each Projective Operator and save the
result as a vector of zeros and ones(though there should only be 1 one). Use the parity matrix to
find the final syndrome, returned as a length 2 vector representing binary. =#
function calculateSyndrome(v:: Array):: Array
    error_vector = []
    push!(error_vector, v' * P_1 * v)
    push!(error_vector, v' * P_2 * v)
    push!(error_vector, v' * P_3 * v)

    return parity_matrix * error_vector
end
#Given a qubit vector v and syndrome vector s, correct qubit vector and return original qubit
function correct(v:: Array, s:: Array):: Array
    syndrome = round(s[1] * 2 + s[2] * 1)
    if syndrome == 0
        return v 
    end
    if syndrome == 1
        return BitFlip1 * v
    end
    if syndrome == 2
        return BitFlip2 * v
    end
    if syndrome == 3
        return BitFlip3 * v
    end
end

#CODE STARTS HERE

#Need to add capabilities for randomizing complex numbers, a little harder and don't want
#to waste too much time on it as it's not super important
alpha_0 = rand()
alpha_1 = sqrt(1 - alpha_0^2)

#Represents state |0> with amplitude alpha_0 and state |1> with amplitude alpha_1 
qubit = [alpha_0, alpha_1]

println("Randomly generated Qubit: ", alpha_0, "|0> + ", alpha_1, "|1>")

#Entangle qubits into state of three qubits
xi = tensor(tensor(qubit, [1 0]), [1 0])

#Encode qubit from |0> -> |000> and |1> -> |111>
encoded_qubit = G_2 * G_1 * xi 
println("Encoded Qubit(Represented as length 8 array): ", encoded_qubit)
#Send encoded qubit through bit flip error channel, new qubit potentially has an error
error_qubit = bitFlipChannel(encoded_qubit, .3)
println("Qubit after bit flip channel: ", error_qubit)

syndrome = calculateSyndrome(error_qubit)
println("Calculated syndrome: ", round(2 * syndrome[1] + syndrome[2]))

corrected_qubit = correct(error_qubit, syndrome)
println("Corrected qubit: ", corrected_qubit)