#ifndef _FUSEDLASSO_ENUMS_H_
#define _FUSEDLASSO_ENUMS_H_

/**
 * @file FusedLasso_enums.h
 * @brief Énumérés centralisés pour les types de régression et de pénalité
 */

/**
 * @enum regEnum
 * @brief Types de régression supportés
 */
enum regEnum {
    GAUSSIAN,  ///< Régression linéaire (Gaussian)
    BINOMIAL   ///< Régression logistique (Binomial)
};

/**
 * @enum penEnum
 * @brief Types de pénalité supportés
 */
enum penEnum {
    L1,        ///< Pénalité L1 (norme L1)
    Huber,     ///< Pénalité de Huber (pénalité robuste)
    L2         ///< Pénalité L2 (norme L2)
};

#endif // _FUSEDLASSO_ENUMS_H_
