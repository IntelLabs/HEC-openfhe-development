//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
CKKS implementation. See https://eprint.iacr.org/2020/1118 for details.
 */

#define PROFILE

#include "cryptocontext.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "scheme/ckksrns/ckksrns-parametergeneration.h"

namespace lbcrypto {

#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
const size_t AUXMODSIZE = 119;
#else
const size_t AUXMODSIZE = 60;
#endif

void ParameterGenerationCKKSRNS::CompositePrimeModuliGen(std::vector<NativeInteger>& moduliQ,
                                                         std::vector<NativeInteger>& rootsQ, usint compositeDegree,
                                                         usint numPrimes, usint firstModSize, usint dcrtBits,
                                                         usint cyclOrder) const {
    std::unordered_set<uint64_t> moduliQRecord;

    // Sample q0, the first primes in the modulus chain
    NativeInteger q;
    uint32_t qBitSize = 0;
    uint32_t remBits  = dcrtBits;
    for (uint32_t d = 1; d <= compositeDegree; ++d) {
        qBitSize = static_cast<uint32_t>(std::ceil(static_cast<double>(remBits) / (compositeDegree - d + 1)));
        q        = FirstPrime<NativeInteger>(qBitSize, cyclOrder);
        q        = PreviousPrime<NativeInteger>(q, cyclOrder);
        while (std::log2(q.ConvertToDouble()) > qBitSize || moduliQRecord.find(q.ConvertToInt()) != moduliQRecord.end())
            q = PreviousPrime<NativeInteger>(q, cyclOrder);
        moduliQ[numPrimes - d] = q;
        rootsQ[numPrimes - d]  = RootOfUnity(cyclOrder, moduliQ[numPrimes - d]);
        moduliQRecord.emplace(q.ConvertToInt());
        remBits -= qBitSize;
    }

#ifdef DEBUG_COMPOSITE_SCALING
    std::cout << __FUNCTION__ << "::" << __LINE__ << " numPrimes=" << numPrimes
              << " compositeDegree=" << compositeDegree << std::endl;
#endif

    std::vector<NativeInteger> qPrev(std::ceil(static_cast<double>(compositeDegree) / 2));
    std::vector<NativeInteger> qNext(compositeDegree - (uint32_t)qPrev.size());

    if (numPrimes > 1) {
        // Prep to compute initial scaling factor
        double sf = moduliQ[numPrimes - 1].ConvertToDouble();
        for (uint32_t d = 2; d <= compositeDegree; ++d) {
            sf *= moduliQ[numPrimes - d].ConvertToDouble();
        }
        double denom = moduliQ[numPrimes - 1].ConvertToDouble();
        for (usint d = 2; d <= compositeDegree; ++d) {
            denom *= moduliQ[numPrimes - d].ConvertToDouble();
        }

        uint32_t cnt = 1;
        for (usint i = numPrimes - compositeDegree; i >= 2 * compositeDegree; i -= compositeDegree) {
            // Compute initial scaling factor
            sf = static_cast<double>(std::pow(sf, 2));
            for (usint d = 0; d < compositeDegree; ++d) {
                sf /= moduliQ[i + d].ConvertToDouble();
            }

            auto sf_sqrt = std::pow((sf * sf) / denom, 1.0 / compositeDegree);  // std::sqrt((sf * sf) / denom);

            NativeInteger sfInt = std::llround(sf_sqrt);
            NativeInteger sfRem = sfInt.Mod(cyclOrder);

#ifdef DEBUG_COMPOSITE_SCALING
            std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " sfInt=" << sfInt
                      << " logq2=" << std::log2(sfInt.ConvertToDouble()) << "\n";
#endif
            double primeProduct = 1.0;
            std::unordered_set<uint64_t> qCurrentRecord;  // current prime tracker

            for (uint32_t step = 0; step < (uint32_t)qPrev.size(); ++step) {
                qPrev[step] = sfInt - (NativeInteger(step + 1) * NativeInteger(cyclOrder)) - sfRem + NativeInteger(1);
                do {
                    qPrev[step] = lbcrypto::PreviousPrime(qPrev[step], cyclOrder);
                } while (moduliQRecord.find(qPrev[step].ConvertToInt()) != moduliQRecord.end() ||
                         qCurrentRecord.find(qPrev[step].ConvertToInt()) != qCurrentRecord.end());
                qCurrentRecord.emplace(qPrev[step].ConvertToInt());
                primeProduct *= qPrev[step].ConvertToDouble();
            }

            for (uint32_t step = 0; step < (uint32_t)qNext.size(); ++step) {
                qNext[step] = sfInt + (NativeInteger(step + 1) * NativeInteger(cyclOrder)) - sfRem + NativeInteger(1);
                do {
                    qNext[step] = lbcrypto::NextPrime(qNext[step], cyclOrder);
                } while (moduliQRecord.find(qNext[step].ConvertToInt()) != moduliQRecord.end() ||
                         qCurrentRecord.find(qNext[step].ConvertToInt()) != qCurrentRecord.end());
                qCurrentRecord.emplace(qNext[step].ConvertToInt());
                primeProduct *= qNext[step].ConvertToDouble();
            }

#ifdef DEBUG_COMPOSITE_SCALING
            std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " primeProduct=" << primeProduct
                      << " sf=" << sf << "\n";
#endif
            if (cnt == 0) {
                NativeInteger qPrevPrev = NativeInteger(qPrev[qPrev.size() - 1].ConvertToInt());
#ifdef DEBUG_COMPOSITE_SCALING
                std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " qPrevPrev=" << qPrevPrev << "\n";
#endif
                while (primeProduct > sf) {
                    do {
                        qCurrentRecord.erase(qPrevPrev.ConvertToInt());  // constant time
                        qPrevPrev = lbcrypto::PreviousPrime(qPrevPrev, cyclOrder);
                    } while (moduliQRecord.find(qPrevPrev.ConvertToInt()) != moduliQRecord.end() ||
                             qCurrentRecord.find(qPrevPrev.ConvertToInt()) != qCurrentRecord.end());
                    qCurrentRecord.emplace(qPrevPrev.ConvertToInt());

#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " primeFound=" << qPrevPrev
                              << " primeProduct=" << primeProduct << " sf=" << sf << "\n";
#endif
                    primeProduct /= qPrev[qPrev.size() - 1].ConvertToDouble();
                    primeProduct *= qPrevPrev.ConvertToDouble();

#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " qPrevPrev" << qPrevPrev
                              << " logq2=" << std::log2(qPrevPrev.ConvertToDouble()) << " primeProduct=" << primeProduct
                              << " sf=" << sf << "\n";
#endif
                }
                for (uint32_t d = 1, p = 0, n = 0; d <= compositeDegree; ++d) {
                    int alternate = d % 2;
                    switch (alternate) {
                        case 1:
                            moduliQ[i - d] = qPrev[p++];
                            break;
                        default:
                            moduliQ[i - d] = qNext[n++];
                            break;
                    }
#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " moduliQ" << moduliQ[i - d]
                              << " logq2=" << std::log2(moduliQ[i - d].ConvertToDouble()) << "\n";
#endif
                    rootsQ[i - d] = RootOfUnity(cyclOrder, moduliQ[i - d]);
                    moduliQRecord.emplace(moduliQ[i - d].ConvertToInt());
                }
                cnt = 1;
            }
            else {
                NativeInteger qNextNext = NativeInteger(qNext[qNext.size() - 1].ConvertToInt());
#ifdef DEBUG_COMPOSITE_SCALING
                std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " qNextNext=" << qNextNext << "\n";
#endif
                while (primeProduct < sf) {
                    do {
                        qCurrentRecord.erase(qNextNext.ConvertToInt());  // constant time
                        qNextNext = lbcrypto::NextPrime(qNextNext, cyclOrder);
                    } while (moduliQRecord.find(qNextNext.ConvertToInt()) != moduliQRecord.end() ||
                             qCurrentRecord.find(qNextNext.ConvertToInt()) != qCurrentRecord.end());
                    qCurrentRecord.emplace(qNextNext.ConvertToInt());

#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " primeFound=" << qNextNext
                              << " primeProduct=" << primeProduct << " sf=" << sf << "\n";
#endif
                    primeProduct /= qNext[qNext.size() - 1].ConvertToDouble();
                    primeProduct *= qNextNext.ConvertToDouble();

#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " qNextNext" << qNextNext
                              << " logq2=" << std::log2(qNextNext.ConvertToDouble()) << " primeProduct=" << primeProduct
                              << " sf=" << sf << "\n";
#endif
                }
                for (uint32_t d = 1, p = 0, n = 0; d <= compositeDegree; ++d) {
                    int alternate = d % 2;
                    switch (alternate) {
                        case 1:
                            moduliQ[i - d] = qPrev[p++];
                            break;
                        default:
                            moduliQ[i - d] = qNext[n++];
                            break;
                    }
#ifdef DEBUG_COMPOSITE_SCALING
                    std::cout << __FUNCTION__ << "::" << __LINE__ << " i=" << i << " moduliQ" << moduliQ[i - d]
                              << " logq2=" << std::log2(moduliQ[i - d].ConvertToDouble()) << "\n";
#endif
                    rootsQ[i - d] = RootOfUnity(cyclOrder, moduliQ[i - d]);
                    moduliQRecord.emplace(moduliQ[i - d].ConvertToInt());
                }
                cnt = 0;
            }
        }  // for loop
    }  // if numPrimes > 1

    if (firstModSize == dcrtBits) {  // this requires dcrtBits < 60
        OPENFHE_THROW(config_error, "firstModSize must be > scalingModSize.");
    }
    else {
        qBitSize = 0;
        remBits  = (uint32_t)firstModSize;
        for (uint32_t d = 1; d <= compositeDegree; ++d) {
            qBitSize = std::ceil(static_cast<double>(remBits) / (compositeDegree - d + 1));
            // Find next prime
            NativeInteger nextInteger = FirstPrime<NativeInteger>(qBitSize, cyclOrder);
            nextInteger               = PreviousPrime<NativeInteger>(nextInteger, cyclOrder);
            // Ensure it fits in 32-bit register
            while (std::log2(nextInteger.ConvertToDouble()) > qBitSize ||
                   moduliQRecord.find(nextInteger.ConvertToInt()) != moduliQRecord.end())
                nextInteger = PreviousPrime<NativeInteger>(nextInteger, cyclOrder);
            // Store prime
            moduliQ[d - 1] = nextInteger;
            rootsQ[d - 1]  = RootOfUnity(cyclOrder, moduliQ[d - 1]);
            // Keep track of existing primes
            moduliQRecord.emplace(moduliQ[d - 1].ConvertToInt());
            remBits -= qBitSize;
#ifdef DEBUG_COMPOSITE_SCALING
            std::cout << __FUNCTION__ << "::" << __LINE__ << " moduliQ=" << moduliQ[d - 1]
                      << " logq2=" << std::log2(moduliQ[d - 1].ConvertToDouble()) << "\n";
#endif
        }
    }

#ifdef DEBUG_COMPOSITE_SCALING
    std::cout << __FUNCTION__ << "::" << __LINE__ << " prime moduli:" << std::endl;
    for (size_t i = 0; i < moduliQ.size(); i++) {
        std::cout << "moduliQ[" << i << "]=" << moduliQ[i] << " logq[" << i
                  << "]=" << std::log2(moduliQ[i].ConvertToDouble()) << std::endl;
    }
#endif

    return;
}

void ParameterGenerationCKKSRNS::SinglePrimeModuliGen(std::vector<NativeInteger>& moduliQ,
                                                      std::vector<NativeInteger>& rootsQ, ScalingTechnique scalTech,
                                                      usint numPrimes, usint firstModSize, usint dcrtBits,
                                                      usint cyclOrder) const {
    NativeInteger q        = FirstPrime<NativeInteger>(dcrtBits, cyclOrder);
    moduliQ[numPrimes - 1] = q;
    rootsQ[numPrimes - 1]  = RootOfUnity(cyclOrder, moduliQ[numPrimes - 1]);

    NativeInteger qNext = q;
    NativeInteger qPrev = q;

    if (numPrimes > 1) {
        if (scalTech != FLEXIBLEAUTO && scalTech != FLEXIBLEAUTOEXT) {
            uint32_t cnt = 0;
            for (usint i = numPrimes - 2; i >= 1; i--) {
                if ((cnt % 2) == 0) {
                    qPrev = lbcrypto::PreviousPrime(qPrev, cyclOrder);
                    q     = qPrev;
                }
                else {
                    qNext = lbcrypto::NextPrime(qNext, cyclOrder);
                    q     = qNext;
                }

                moduliQ[i] = q;
                rootsQ[i]  = RootOfUnity(cyclOrder, moduliQ[i]);
                cnt++;
            }
        }
        else {
            double sf    = moduliQ[numPrimes - 1].ConvertToDouble();
            uint32_t cnt = 1;
            for (usint i = numPrimes - 2; i >= 1; i--) {
                sf = static_cast<double>(pow(sf, 2) / moduliQ[i + 1].ConvertToDouble());
                if ((cnt % 2) == 0) {
                    NativeInteger sfInt = std::llround(sf);
                    NativeInteger sfRem = sfInt.Mod(cyclOrder);
                    NativeInteger qPrev = sfInt - NativeInteger(cyclOrder) - sfRem + NativeInteger(1);

                    bool hasSameMod = true;
                    while (hasSameMod) {
                        hasSameMod = false;
                        qPrev      = lbcrypto::PreviousPrime(qPrev, cyclOrder);
                        // qPrev      = lbcrypto::PreviousPrime(qPrev - 1000000000000000, cyclOrder);
                        for (uint32_t j = i + 1; j < numPrimes; j++) {
                            if (qPrev == moduliQ[j]) {
                                hasSameMod = true;
                            }
                        }
                    }
                    moduliQ[i] = qPrev;
                }
                else {
                    NativeInteger sfInt = std::llround(sf);
                    NativeInteger sfRem = sfInt.Mod(cyclOrder);
                    NativeInteger qNext = sfInt + NativeInteger(cyclOrder) - sfRem + NativeInteger(1);

                    bool hasSameMod = true;
                    while (hasSameMod) {
                        hasSameMod = false;
                        qNext      = lbcrypto::NextPrime(qNext, cyclOrder);
                        for (uint32_t j = i + 1; j < numPrimes; j++) {
                            if (qNext == moduliQ[j]) {
                                hasSameMod = true;
                            }
                        }
                    }
                    moduliQ[i] = qNext;
                }

                rootsQ[i] = RootOfUnity(cyclOrder, moduliQ[i]);
                cnt++;
            }
        }
    }

    if (firstModSize == dcrtBits) {  // this requires dcrtBits < 60
        OPENFHE_THROW(config_error, "firstModSize must be > scalingModSize.");
        moduliQ[0] = PreviousPrime<NativeInteger>(qPrev, cyclOrder);
    }
    else {
        NativeInteger firstInteger = FirstPrime<NativeInteger>(firstModSize, cyclOrder);
        moduliQ[0]                 = PreviousPrime<NativeInteger>(firstInteger, cyclOrder);
    }

    rootsQ[0] = RootOfUnity(cyclOrder, moduliQ[0]);
}

bool ParameterGenerationCKKSRNS::ParamsGenCKKSRNS(std::shared_ptr<CryptoParametersBase<DCRTPoly>> cryptoParams,
                                                  usint cyclOrder, usint numPrimes, usint scalingModSize,
                                                  usint firstModSize, uint32_t numPartQ,
                                                  COMPRESSION_LEVEL mPIntBootCiphertextCompressionLevel) const {
    const auto cryptoParamsCKKSRNS = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cryptoParams);

    KeySwitchTechnique ksTech        = cryptoParamsCKKSRNS->GetKeySwitchTechnique();
    ScalingTechnique scalTech        = cryptoParamsCKKSRNS->GetScalingTechnique();
    EncryptionTechnique encTech      = cryptoParamsCKKSRNS->GetEncryptionTechnique();
    MultiplicationTechnique multTech = cryptoParamsCKKSRNS->GetMultiplicationTechnique();
    ProxyReEncryptionMode PREMode    = cryptoParamsCKKSRNS->GetPREMode();

    // Determine appropriate composite degree automatically if scaling technique set to COMPOSITESCALINGAUTO
    cryptoParamsCKKSRNS->ConfigureCompositeDegree(firstModSize);
    usint compositeDegree = cryptoParamsCKKSRNS->GetCompositeDegree();
    // compositeDegree *= (usint)1;  // @fdiasmor: Avoid unused variable compilation error.
    // if (compositeDegree > 2 && firstModSize > 64) {
    //     OPENFHE_THROW(config_error, "This COMPOSITESCALING* version does not support composite degree > 2.");
    // }
    // Bookeeping unique prime moduli
    //    std::unordered_set<uint64_t> moduliQRecord;

    if ((PREMode != INDCPA) && (PREMode != NOT_SET)) {
        std::stringstream s;
        s << "This PRE mode " << PREMode << " is not supported for CKKSRNS";
        OPENFHE_THROW(s.str());
    }

    usint extraModSize = 0;
    if (scalTech == FLEXIBLEAUTOEXT) {
        // TODO: Allow the user to specify this?
        extraModSize = DCRT_MODULUS::DEFAULT_EXTRA_MOD_SIZE;
    }

    //// HE Standards compliance logic/check
    SecurityLevel stdLevel = cryptoParamsCKKSRNS->GetStdLevel();
    uint32_t auxBits       = (scalTech == COMPOSITESCALINGAUTO || scalTech == COMPOSITESCALINGMANUAL) ? 30 : AUXMODSIZE;
    uint32_t n             = cyclOrder / 2;
    uint32_t qBound        = firstModSize + (numPrimes - 1) * scalingModSize + extraModSize;
    // Estimate ciphertext modulus Q bound (in case of GHS/HYBRID P*Q)
    if (ksTech == HYBRID) {
        if (scalTech == COMPOSITESCALINGAUTO || scalTech == COMPOSITESCALINGMANUAL) {
            uint32_t tmpFactor = (compositeDegree == 2) ? 2 : 4;
            qBound += ceil(ceil(static_cast<double>(qBound) / numPartQ) / (tmpFactor * auxBits)) * tmpFactor * auxBits;
        }
        else {
            qBound += ceil(ceil(static_cast<double>(qBound) / numPartQ) / auxBits) * auxBits;
        }
    }

    // GAUSSIAN security constraint
    DistributionType distType = (cryptoParamsCKKSRNS->GetSecretKeyDist() == GAUSSIAN) ? HEStd_error : HEStd_ternary;
    auto nRLWE                = [&](usint q) -> uint32_t {
        return StdLatticeParm::FindRingDim(distType, stdLevel, q);
    };

    // Case 1: SecurityLevel specified as HEStd_NotSet -> Do nothing
    if (stdLevel != HEStd_NotSet) {
        if (n == 0) {
            // Case 2: SecurityLevel specified, but ring dimension not specified

            // Choose ring dimension based on security standards
            n         = nRLWE(qBound);
            cyclOrder = 2 * n;
        }
        else {  // if (n!=0)
            // Case 3: Both SecurityLevel and ring dimension specified

            // Check whether particular selection is standards-compliant
            auto he_std_n = nRLWE(qBound);
            if (he_std_n > n) {
                OPENFHE_THROW("The specified ring dimension (" + std::to_string(n) +
                              ") does not comply with HE standards recommendation (" + std::to_string(he_std_n) + ").");
            }
        }
    }
    else if (n == 0) {
        OPENFHE_THROW("Please specify the ring dimension or desired security level.");
    }
    //// End HE Standards compliance logic/check

    usint dcrtBits = scalingModSize;

    if (scalTech == COMPOSITESCALINGAUTO || scalTech == COMPOSITESCALINGMANUAL) {
        numPrimes *= compositeDegree;
    }

    uint32_t vecSize = (extraModSize == 0) ? numPrimes : numPrimes + 1;
    std::vector<NativeInteger> moduliQ(vecSize);
    std::vector<NativeInteger> rootsQ(vecSize);

    if (scalTech == COMPOSITESCALINGAUTO || scalTech == COMPOSITESCALINGMANUAL) {
        CompositePrimeModuliGen(moduliQ, rootsQ, compositeDegree, numPrimes, firstModSize, dcrtBits, cyclOrder);
    }
    else {
        SinglePrimeModuliGen(moduliQ, rootsQ, scalTech, numPrimes, firstModSize, dcrtBits, cyclOrder);
    }

    if (scalTech == FLEXIBLEAUTOEXT) {
        // no need for extra checking as extraModSize is automatically chosen by the library
        moduliQ[numPrimes] = FirstPrime<NativeInteger>(extraModSize - 1, cyclOrder);
        rootsQ[numPrimes]  = RootOfUnity(cyclOrder, moduliQ[numPrimes]);
    }

    auto paramsDCRT = std::make_shared<ILDCRTParams<BigInteger>>(cyclOrder, moduliQ, rootsQ);

    cryptoParamsCKKSRNS->SetElementParams(paramsDCRT);

    const EncodingParams encodingParams = cryptoParamsCKKSRNS->GetEncodingParams();
    if (encodingParams->GetBatchSize() > n / 2)
        OPENFHE_THROW("The batch size cannot be larger than ring dimension / 2.");

    if (encodingParams->GetBatchSize() & (encodingParams->GetBatchSize() - 1))
        OPENFHE_THROW("The batch size can only be set to zero (for full packing) or a power of two.");

    // if no batch size was specified, we set batchSize = n/2 by default (for full
    // packing)
    if (encodingParams->GetBatchSize() == 0) {
        uint32_t batchSize = n / 2;
        EncodingParams encodingParamsNew(
            std::make_shared<EncodingParamsImpl>(encodingParams->GetPlaintextModulus(), batchSize));
        cryptoParamsCKKSRNS->SetEncodingParams(encodingParamsNew);
    }

    cryptoParamsCKKSRNS->PrecomputeCRTTables(ksTech, scalTech, encTech, multTech, numPartQ, auxBits, extraModSize);

    return true;
}

}  // namespace lbcrypto
