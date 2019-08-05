#include <regex>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/find_iterator.hpp>
#include "Phreeqc.h"
#include "phqalloc.h"


/* ---------------------------------------------------------------------- */
void Phreeqc::parse_eq(std::string reaction_equation, struct elt_list** elt_ptr,
                       int association)
/* ---------------------------------------------------------------------- */
/*
 *   function to break equation up into component species
 *   returns species name, coefficient, and charge in global variable "rxn_list".
 *   Also returns pointer to elt_list structure array with list of elements
 *   and coefficients in species.
 *   rxn_list     global variable, output with contiguous rxn_list structures,
 *                      which is the reaction.  Target species is first in
 *                      list, coefficient is -1.0.
 *
 *   Argurments:
 *      *eqn         input,  pointer to equation to be parsed.
 *     **elt_ptr     output, pointer to contiguous elt_list structures,
 *                           which is list of elements and coefficients in the
 *                           target species.
 *       association input, TRUE or FALSE if reaction is an association reaction
 */
{
    int i;
    LDBLE coef, l_z;
    char c;
    std::string ptr;
    std::string token;

    paren_count = 0;
    // split reaction equation into two parts
    std::string delimiter = "=";
    std::string reaction_equation_lhs =
        reaction_equation.substr(0, reaction_equation.find(delimiter));
    std::string reaction_equation_rhs =
        reaction_equation.substr(reaction_equation.find(delimiter) + 1);

    // Parse stoichiometric coefficients, species names, and electric charges
    // carried
    parseReactionEquationlhs(reaction_equation_lhs);
    parseReactionEquationrhs(reaction_equation_rhs);

    /*
     *   Sort list of reaction species
     */
//    trxn_sort();
    /*
     *   Get elements in species or mineral formula
     */
    //	count_elts = 0;
    //	replace("(s)", "", token);
    //	replace("(S)", "", token);
    //	replace("(g)", "", token);
    //	replace("(G)", "", token);

    /*
     *   Sort elements in reaction and combine
     */
    //	qsort(elt_list, (size_t) count_elts, (size_t) sizeof(struct elt_list),
    //		  elt_list_compare);
    //	if (elt_list_combine() == ERROR)
    //		return (ERROR);
    /*
     *   Malloc space and store element data for return
     */
//        *elt_ptr =
//            (struct elt_list *) PHRQ_malloc((size_t) (count_elts + 1) *
//                                            sizeof(struct elt_list));
    //	if (*elt_ptr == NULL)
    //	{
    //		malloc_error();
    //	}
    //	else
    //	{
    //		for (i = 0; i < count_elts; i++)
    //		{
    //			(*elt_ptr)[i].elt = elt_list[i].elt;
    //			(*elt_ptr)[i].coef = -elt_list[i].coef;
    //		}
    //		(*elt_ptr)[count_elts].elt = NULL;
    //	}
    /*
     *   Debugging print of parsed equation
        trxn_print();
     */
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_eqn(int association)
/* ---------------------------------------------------------------------- */
/*
 *   Check that equation is balanced reaction. Uses array "trxn.token" and
 *   assumes "count_trxn" is set to number of species in reaction.
 *   Charge and elements are checked.
 *
 *   Arguments:
 *
 *      association  input, TRUE or FALSE if reaction is an association reaction.
 */
{
    int i;
    int oops = 0;
    LDBLE sumcharge;

    paren_count = 0;
    count_elts = 0;
/*
 *   Check that coefficient of first species is -1.0
 */
    if (equal(trxn.token[0].coef, -1.0, TOL) == FALSE)
    {
        if (association == TRUE)
        {
            error_string = sformatf(
                    "Coefficient of first species on rhs is not equal to 1.0.");
            error_msg(error_string, CONTINUE);
        }
        else
        {
            error_string = sformatf(
                    "Coefficient of mineral (first on lhs) is not equal to 1.0.");
            error_msg(error_string, CONTINUE);
        }
        return (ERROR);
    }
/*
 *   Go through all species in the reaction; sum the charge and store elements
 */
    sumcharge = 0.0;
    for (i = 0; i < count_trxn; i++)
    {
        sumcharge += (trxn.token[i].coef) * (trxn.token[i].charge);
        std::string temp_name = trxn.token[i].name;
        std::string t_ptr = temp_name;
        if (get_elts_in_species(t_ptr, trxn.token[i].coef) == ERROR)
        {
//			free_check_null(temp_name);
            return (ERROR);
        }
//		free_check_null(temp_name);
    }
/*
 *   Sort elements in reaction and combine
 */
    //    qsort(elt_list, (size_t) count_elts, (size_t) sizeof(struct elt_list),
    //          elt_list_compare);
    if (elt_list_combine() == ERROR)
        return (ERROR);
/*
 *   Check charge
 */
    if (equal(sumcharge, 0.0, TOL) == FALSE)
    {
        error_string = sformatf( "Equation is not charge balanced, right - left = %7.4f moles charge", sumcharge);
        error_msg(error_string, CONTINUE);
        oops++;
    }
/*
 *   Check mass balance
 */
    for (i = 0; i < count_elts; i++)
    {
        if ((equal(elt_list[i].coef, 0.0, TOL) == FALSE) /*&&
            strncmp((elt_list[i].elt).name, "e", MAX_LENGTH) != 0*/)
        {
            error_string = sformatf(
                "Equation does not balance for element, %s: right - left = %7.4f moles",
                (elt_list[i].elt).name, elt_list[i].coef);
            error_msg(error_string, CONTINUE);
            oops++;
        }
    }
    if (oops == 0)
    {
        return (OK);
    }
    else
    {
        return (ERROR);
    }
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_charge(char *charge, LDBLE * l_z)
/* ---------------------------------------------------------------------- */
/*
 *   Function takes character string and calculates the charge on
 *   the species.  Charge can be in two forms: (1) string of "+" or "-"
 *   or (2) + or - followed by an integer. Charge is reduced to form (2)
 *   and stored in the pointer location *charge.
 *
 *   Arguments:
 *      *charge    input, string containing charge
 *                 output, string containing charge of the form + or -
 *                    followed by an integer, if integer greater than 1.
 *      *z         output, value of charge.
 *
 *   Returns:
 *      ERROR,
 *      OK.
 */
{
    int i;
    char *ptr;
    char c, c1;
/*
 *   Charge is zero
 */
    if ((c = charge[0]) == '\0')
    {
        *l_z = 0.0;
        return (OK);
    }
/*
 *   Error check for + or - at start of string
 */
    if (c != '+' && c != '-')
    {
        error_string = sformatf(
                "Character string for charge does not start with + or -,\t%s.",
                charge);
        error_msg(error_string, CONTINUE);
        return (ERROR);
    }
/*
 *   Count string of +'s or -'s
 */
    i = 0;
    while (c == (c1 = charge[i++]));
    i--;
    if (c1 == '\0')
    {
        if (c == '-')
            i = -i;
    }
    else
    {
/*
 *   + or - followed by a number
 */
        errno = 0;
        i = strtol(charge, &ptr, 0);
/*
 *   Truncate fractional part of charge if all zeros
 */
        if (*ptr != '\0')
        {
            if (*ptr == '.')
            {
                while (*(++ptr) != '\0')
                {
                    if (*ptr != '0')
                    {
                        *l_z = strtod(charge, &ptr);
                        return (OK);
                    }
                }
/*
 *   Non-numeric characters
 */
            }
            else
            {
                error_string = sformatf(
                        "Error in character string for charge, %s.", charge);
                error_msg(error_string, CONTINUE);
                return (ERROR);
            }
        }
    }
/*
 *   Charge is zero, must have had +0 or -0 in eqn
 */
    if (i == 0)
    {
        charge[0] = '\0';
    }
/*
 *   Charge is +1 or -1, single + or -
 */
    if (abs(i) == 1)
    {
        charge[0] = c;
        charge[1] = '\0';
    }
/*
 *   Abs(z)>1, set charge to + or - plus integer
 */
    if (abs(i) > 1)
    {
        if (sprintf(charge, "%-+d", i) == EOF)
        {
            error_string = sformatf(
                    "Error converting charge to character string, %s.",
                    charge);
            error_msg(error_string, CONTINUE);
            return (ERROR);
        }
    }
    *l_z = i;
    return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_coef(LDBLE * coef, std::string eqnaddr)
/* ---------------------------------------------------------------------- */
/*
 *   Function reads through eqn and determines the coefficient of the next
 *   species.
 *
 *   Arguments:
 *       *coef    output, coefficient of next species.
 *
 *      **eqnaddr input, pointer to a position in eqn to start parsing
 *                output, pointer to next position after coefficient
 *
 *   Returns:
 *      ERROR,
 *      OK.
 */
{
    int i;
    char c, c1;
    char *ptr, *ptr1, *rest;
    std::string token;

//	rest = *eqnaddr;
//	ptr = *eqnaddr;				/* address of a position in eqn */
    c = *ptr;					/* character in eqn */
    *coef = 0.0;
/*
 *   No leading sign or number
 */
    if (isalpha((int) c) ||
        (c == '(') || (c == ')') || (c == '[') || (c == ']'))
    {
        *coef = 1.0;
        return (OK);
    }
/*
 *   Leading +, no digits
 */
    c1 = *(ptr + 1);
    if (c == '+' &&
        (isalpha((int) c1) ||
         (c1 == '(') || (c1 == ')') || (c1 == '[') || (c1 == ']')))
    {
//		*eqnaddr = ++ptr;
        *coef = 1.0;
        return (OK);
    }
/*
 *   Leading -, no digits
 */
    if (c == '-' &&
        (isalpha((int) c1) ||
         (c1 == '(') || (c1 == ')') || (c1 == '[') || (c1 == ']')))
    {
//		*eqnaddr = ++ptr;
        *coef = -1.0;
        return (OK);
    }
    i = 0;
/*
 *   Has number coefficient
 */
    if (isdigit((int) c) || c == '+' || c == '-' || c == '.')
    {
        while (isdigit((int) c) || c == '+' || c == '-' || c == '.')
        {
            token[i++] = c;
//			if (i >= MAX_LENGTH)
//			{
//				error_string = sformatf(
//						"Coefficient has more than MAX_LENGTH characters.");
//				error_msg(error_string, CONTINUE);
//				return (ERROR);
//			}
            c = *(++ptr);
        }
        token[i] = '\0';
//		*eqnaddr = ptr;
        errno = 0;
        *coef = std::stod(token);
        if ((errno == ERANGE) || (*ptr1 != '\0'))
        {
            error_string = sformatf(
                    "Error converting coefficient in get_coef, %s.", token);
            error_msg(error_string, CONTINUE);
            return (ERROR);
        }
        return (OK);
    }
/*
 *   None of the above, unknown construct
 */
    error_string = sformatf(
            "Illegal equation construct detected in get_coef.\n\t%s.", rest);
    error_msg(error_string, CONTINUE);
    return (ERROR);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_elt(std::string t_ptr, std::string element, int *i)
/* ---------------------------------------------------------------------- */
/*
 *      Function reads an element name out of the equation string.
 *      An element name is composed of a capital letter followed by any number
 *      of lower case characters.
 *
 *      Arguments:
 *         **t_ptr   input, points to position in the equation to begin
 *                   output, points to next character of equation after
 *                      element name.
 *         *element  input pointer to place to return element character string
 */
{
    char c;

//	c = *(*t_ptr)++;
//	if (c == '\0')
//	{
//		error_string = sformatf(
//				"Empty string in get_elt.  Expected an element name.");
//		error_msg(error_string, CONTINUE);
//		return (ERROR);
//	}
/*
// *   Load name into char array element
// */
//	element[0] = c;
//	*i = 1;
//	if (c == '[')
//	{
//		while ((c = (**t_ptr)) != ']')
//		{
//			element[*i] = c;
//			(*i)++;
//			(*t_ptr)++;
//			if ((c = (**t_ptr)) == ']')
//			{
//				element[*i] = c;
//				(*i)++;
//				(*t_ptr)++;
//				break;
//			}
//			else if (**t_ptr == '\0')
//			{
//				error_msg("No ending bracket (]) for element name", CONTINUE);
//				input_error++;
//				break;
//			}
//		}
//		while (islower((int) (c = (**t_ptr))) || c == '_')
//		{
//			element[*i] = c;
//			(*i)++;
//			(*t_ptr)++;
//		}
//	}
//	else
//	{
//		while (islower((int) (c = (**t_ptr))) || c == '_')
//		{
//			element[*i] = c;
//			(*i)++;
//			(*t_ptr)++;
//		}
//	}
    element[*i] = '\0';
    return (OK);
}

std::vector<std::string> s_split(std::string const& str)
{
    std::string str1 = str;
    std::string::size_type identifier_pos;
    if ((identifier_pos = str.find_first_of("+-")) != std::string::npos)
        str1 = str.substr(0, identifier_pos);

    std::vector<std::string> match_str;
    std::string const delims{ "ABCDEFGHIJKLMNOPQRSTUVWXYZ" };

    size_t beg, pos = 0;
    while ((beg = str1.find_first_of(delims, pos)) != std::string::npos)
    {
        pos = str1.find_first_of(delims, beg + 1);
        std::string element_name = str1.substr(beg, pos - beg);
//        elt_list.emplace_back(element_name, coefficient);
        match_str.push_back(str1.substr(beg, pos - beg));
    }
    return match_str;
}

int Phreeqc::
get_elts_in_species(std::string species, double coef)
/* ---------------------------------------------------------------------- */
{
/*
 *    Makes a list of elements with their coefficients, stores elements
 *    in elt_list at position count_elts.  Global variable count_elts is
 *    updated with each stored element.  Also uses static global variable
 *    paren_count.
 *
 *    Arguments:
 *       **t_ptr    input, point in token string to start looking
 *                  output, is next position to start looking
 *         coef     input, coefficient to multiply subscripts by
 */
    double d;
    std::vector<std::string> elements;
    if (species == "e-")
    {
        elements.push_back("e");
        d = 1.0;
    }
    else
    {
        elements = s_split(species);
    }

    //	while (((c = **t_ptr) != '+') && (c != '-') && (c != '\0'))
//    auto coefficeint = d * coef;

    //			if (get_num(t_ptr, &d) == ERROR)
    //			{
    //				return (ERROR);
    //			}
    //			count_elts++;
    /*
    // *   Expand working space for elements if necessary
    // */
    //			if (count_elts >= max_elts)
    //			{
    //				space((void **) ((void *) &elt_list), count_elts, &max_elts,
    //					  sizeof(struct elt_list));
    //			}
    //			continue;
    //		}
    /*
    // *   Open parentheses
    // */
    //		if (c == '(')
    //		{
    //			count = count_elts;
    //			if (c1 == ')')
    //			{
    //				error_string = sformatf( "Empty parentheses.");
    //				warning_msg(error_string);
    //			}
    //			paren_count++;
    //			(*t_ptr)++;
    //			if (get_elts_in_species(t_ptr, coef) == ERROR)
    //			{
    //				return (ERROR);
    //			}
    //			if (get_num(t_ptr, &d) == ERROR)
    //			{
    //				return (ERROR);
    //			}
    //			for (i = count; i < count_elts; i++)
    //			{
    //				elt_list[i].coef *= d;
    //			}
    //			continue;
    //		}
    /*
    // *   Colon
    // */
    //		if (c == ':')
    //		{
    //			count = count_elts;
    //			(*t_ptr)++;
    //			if (get_num(t_ptr, &d) == ERROR)
    //			{
    //				return (ERROR);
    //			}
    //			if (get_elts_in_species(t_ptr, coef) == ERROR)
    //			{
    //				return (ERROR);
    //			}
    //			for (i = count; i < count_elts; i++)
    //			{
    //				elt_list[i].coef *= d;
    //			}
    //			continue;
    //		}
    /*
    // *   Not beginning of element and not opening paren
    // */
    //		error_string = sformatf(
    //				"Parsing error in get_elts_in_species, unexpected character,
    //%c.", 				c); 		error_msg(error_string, CONTINUE); 		input_error++; 		return
    //(ERROR);
    //	}
    if (paren_count != 0)
    {
        error_string = sformatf( "Unbalanced parentheses.");
        error_msg(error_string, CONTINUE);
        input_error++;
        return (ERROR);
    }
    return (OK);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
get_secondary(char **t_ptr, char *element, int *i)
/* ---------------------------------------------------------------------- */
/*
 *      Function reads an element name out of the equation string.
 *      An element name is composed of a capital letter followed by any number
 *      of lower case characters.
 *
 *      Arguments:
 *         **t_ptr   input, points to position in the equation to begin
 *                   output, points to next character of equation after
 *                      element name.
 *         *element  input pointer to place to return element character string
 */
{
    int j;
    char c;
    char *ptr;

    c = *(*t_ptr)++;
    if (c == '\0')
    {
        error_string = sformatf(
                "Empty string in get_elt.  Expected an element name.");
        error_msg(error_string, CONTINUE);
        input_error++;
        return (ERROR);
    }
/*
 *   Load name into char array element
 */
    element[0] = c;
    *i = 1;
    if (c == '[')
    {
        while ((c = (**t_ptr)) != ']')
        {
            element[*i] = c;
            (*i)++;
            (*t_ptr)++;
            if ((c = (**t_ptr)) == ']')
            {
                element[*i] = c;
                (*i)++;
                (*t_ptr)++;
                c = (**t_ptr);
                break;
            }
            else if ((c = (**t_ptr)) == '\0')
            {
                error_msg("Did not find ending bracket (])", CONTINUE);
                input_error++;
                return (ERROR);
            }
        }
        while (islower((int) (c = (**t_ptr))) || c == '_')
        {
            element[*i] = c;
            (*i)++;
            (*t_ptr)++;
        }
    }
    else
    {
        while (islower((int) (c = (**t_ptr))) || c == '_')
        {
            element[*i] = c;
            (*i)++;
            (*t_ptr)++;
        }
    }
/*
 *   Check if secondary master species element
 */
    j = *i;
    ptr = *t_ptr;
    if (c == '(')
    {
        /* copy parenthesis */
        element[*i] = c;
        (*i)++;
        (*t_ptr)++;
        /* copy number */
        for (;;)
        {
            c = **t_ptr;
            if (isdigit((int) c) || c == '-' || c == '.')
            {
                element[*i] = c;
                (*i)++;
                (*t_ptr)++;
            }
            else if (c == '+')
            {
                (*t_ptr)++;
            }
            else
            {
                break;
            }
        }
        /* go back to before parenthesis */
        if (c != ')')
        {
            *i = j;
            *t_ptr = ptr;
            /* put in closing parenthesis */
        }
        else
        {
            element[*i] = c;
            (*i)++;
            (*t_ptr)++;
        }
    }
    element[*i] = '\0';
    return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_secondary_in_species(std::string t_ptr, LDBLE coef)
/* ---------------------------------------------------------------------- */
{
/*
 *    Makes a list of elements with their coefficients, stores elements
 *    in elt_list at position count_elts.  Global variable count_elts is
 *    updated with each stored element.  Also uses static global variable
 *    paren_count.
 *
 *    Arguments:
 *       **t_ptr    input, point in token string to start looking
 *                  output, is next position to start looking
 *         coef     input, coefficient to multiply subscripts by
 */
    int i, count, l;
    char c, c1;
    LDBLE d;
    std::string element;

//	while (((c = **t_ptr) != '+') && (c != '-') && (c != '\0'))
//	{
//		/* close parenthesis */
//		if (c == ')')
//		{
//			paren_count--;
//			if (paren_count < 0)
//			{
//				error_string = sformatf( "Too many right parentheses.");
//				error_msg(error_string, CONTINUE);
//				input_error++;
//				return (ERROR);
//			}
//			(*t_ptr)++;
//			return (OK);
//		}
//		c1 = *((*t_ptr) + 1);
//		/* beginning of element name */
//		if (isupper((int) c) || c == '[' || (c == 'e' && c1 == '-'))
//		{
///*
// *   Get new element and subscript
// */
//			if (get_secondary(t_ptr, element, &l) == ERROR)
//			{
//				return (ERROR);
//			}
//			elt_list[count_elts].elt = element_store(element);
//			if (get_num(t_ptr, &d) == ERROR)
//			{
//				return (ERROR);
//			}
//			elt_list[count_elts].coef = d * coef;
//			count_elts++;
///*
// *   Expand working space for elements if necessary
// */
//			if (count_elts >= max_elts)
//			{
//				space((void **) ((void *) &elt_list), count_elts, &max_elts,
//					  sizeof(struct elt_list));
//			}
//			continue;
//		}
///*
// *   Open parentheses
// */
//		if (c == '(')
//		{
//			count = count_elts;
//			if (c1 == ')')
//			{
//				error_string = sformatf( "Empty parentheses.");
//				warning_msg(error_string);
//			}
//			paren_count++;
//			(*t_ptr)++;
//			if (get_secondary_in_species(t_ptr, coef) == ERROR)
//			{
//				return (ERROR);
//			}
//			if (get_num(t_ptr, &d) == ERROR)
//			{
//				return (ERROR);
//			}
//			for (i = count; i < count_elts; i++)
//			{
//				elt_list[i].coef *= d;
//			}
//			continue;
//		}
///*
// *   Colon
// */
//		if (c == ':')
//		{
//			count = count_elts;
//			(*t_ptr)++;
//			if (get_num(t_ptr, &d) == ERROR)
//			{
//				return (ERROR);
//			}
//			if (get_secondary_in_species(t_ptr, coef) == ERROR)
//			{
//				return (ERROR);
//			}
//			for (i = count; i < count_elts; i++)
//			{
//				elt_list[i].coef *= d;
//			}
//			continue;
//		}
///*
// *   Not beginning of element and not opening paren
// */
//		error_string = sformatf(
//				"Parsing error in get_secondary_in_species, unexpected character, %c.",
//				c);
//		error_msg(error_string, CONTINUE);
//		return (ERROR);
//	}
    if (paren_count != 0)
    {
        error_string = sformatf( "Unbalanced parentheses.");
        error_msg(error_string, CONTINUE);
        return (ERROR);
    }
    return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_num(std::string t_ptr, LDBLE * num)
/* ---------------------------------------------------------------------- */
/*
 *      Function reads through a string looking for leading numeric field
 *      if no numeric field is found, the number is set to 1.0.
 *
 *      Arguments:
 *
 *       **t_ptr     input, points to a position in a character string from which
 *                      a number is to be extracted.
 *                   output, points to next position to be parsed.
 *        *num       address where the number is to be stored.
 *
 *      Returns:
 *         ERROR,
 *         OK.
 */
{
    int i, decimal;
    char c;
    char *ptr1;
    std::string token;

    *num = 1.0;
    i = 0;
//	c = **t_ptr;
    decimal = 0;
    if (isdigit((int) c) || (c == '.'))
    {
        while (isdigit((int) c) || (c == '.'))
        {
            if (c == '.')
                decimal++;
            if (decimal > 1)
                break;
            token[i++] = c;
            /* check number length */
//			if (i >= MAX_LENGTH)
//			{
//				error_string = sformatf(
//						"Number was greater than MAX_LENGTH characters.");
//				error_msg(error_string, CONTINUE);
//				input_error++;
//				return (ERROR);
//			}
//			c = *(++(*t_ptr));
        }
        token[i] = '\0';
        errno = 0;
        *num = std::stod(token);
        if (errno == ERANGE)
        {
            error_string = sformatf( "Converting number in get_num, %s.", token);
            input_error++;
            error_msg(error_string, CONTINUE);
            return (ERROR);
        }
    }
    return (OK);
}

std::pair<std::string, double> Phreeqc::parseSpeciesName(std::string term)
{
    std::array<std::string, 2> identifiers = {"+", "-"};
    int charge = 0;
    int occurrences = 0;
    for (auto const& identifier : identifiers)
    {
        std::string::size_type start = 0;
        while ((start = term.find(identifier, start)) != std::string::npos)
        {
            ++occurrences;
            start += identifier.length();
        }

        if (occurrences == 0)
            continue;
        else if (occurrences == 1 && term.back() == identifier[0])
        {
            assert(term.find_last_of(identifier) -
                       term.find_first_of(identifier) + 1 ==
                   occurrences);
            charge = identifier == "+" ? 1 : -1;
            return std::make_pair(term, charge);
        }
        else if (occurrences == 1 && std::isdigit(term.back()))
        {
            assert(term.find_last_of(identifier) -
                       term.find_first_of(identifier) + 1 ==
                   occurrences);
            charge = static_cast<int>(term.back()) - '0';
            charge *= identifier == "+" ? 1 : -1;
            return std::make_pair(term, charge);
        }
        else
        {
            assert(term.find_last_of(identifier) -
                       term.find_first_of(identifier) + 1 ==
                   occurrences);
            charge = occurrences;
            charge *= identifier == "+" ? 1 : -1;
            return std::make_pair(term.substr(0, term.find(identifier) + 1) +
                      std::to_string(occurrences), charge);
        }
    }

    if (occurrences == 0)
    {
        return std::make_pair(term, 0);
    }
}

struct species& Phreeqc::getOrcreateSpecies(std::pair<std::string, double> species_name)
{
    auto& name = species_name.first;
    auto& charge = species_name.second;
    auto compare_by_name = [&name](auto const& species) {
        return species.name == name;};
    auto it = std::find_if(species_set.begin(), species_set.end(), compare_by_name);
    if (it == species_set.end())
    {
        species_set.emplace_back(name, charge);
        return species_set.back();
    }
    else
    {
        return *it;
    }
}

void Phreeqc::parseReactionEquationlhs(std::string& expression)
{
    // split expression into terms
    std::vector<std::string> terms;
    typedef boost::split_iterator<std::string::iterator> string_split_iterator;
    for (string_split_iterator It = boost::make_split_iterator(
             expression, boost::first_finder(" +", boost::is_iequal()));
         It != string_split_iterator();
         ++It)
    {
        auto term = boost::copy_range<std::string>(*It);
        // compress the terms trimmed by removing whitespace
        boost::erase_all(term, " ");
        terms.push_back(term);
    }

    for (auto& term : terms)
    {
        // extract stoichiometric coefficient
        std::string::size_type nondigit_position =
            term.find_first_not_of("0123456789.");
        double coefficient = nondigit_position == 0
                                 ? 1.0
                                 : std::stod(term.substr(0, nondigit_position));
        term.erase(0, nondigit_position);

        // parse species names and electronic charges carried
        auto species = parseSpeciesName(term);
        auto& species_name = species.first;
        auto& charge = species.second;

        // initialization
        trxn.token.emplace_back(species_name, charge, coefficient);
    }
}

void Phreeqc::parseReactionEquationrhs(std::string& expression)
{
    // split expression into terms
    std::vector<std::string> terms;
    typedef boost::split_iterator<std::string::iterator> string_split_iterator;
    for (string_split_iterator It = boost::make_split_iterator(
             expression, boost::first_finder(" +", boost::is_iequal()));
         It != string_split_iterator();
         ++It)
    {
        auto term = boost::copy_range<std::string>(*It);
        // compress the terms trimmed by removing whitespace
        boost::erase_all(term, " ");
        terms.push_back(term);
    }

    int i = 0;
    for (auto& term : terms)
    {
        // extract stoichiometric coefficient
        std::string::size_type nondigit_position =
            term.find_first_not_of("0123456789.");
        double coefficient = nondigit_position == 0
                                 ? -1.0
                                 : -1.0 * std::stod(term.substr(0, nondigit_position));
        term.erase(0, nondigit_position);

        // parse species names and electronic charges carried
        auto species = parseSpeciesName(term);
        auto& species_name = species.first;
        auto& charge = species.second;

        // initialization
        if (i == 0)
        {
            auto it = trxn.token.begin();
            trxn.token.emplace(it, species_name, charge, coefficient);
            get_elts_in_species(species_name, coefficient);
        }
        else
        {
            trxn.token.emplace_back(species_name, charge, coefficient);
        }
        ++i;
    }
}
