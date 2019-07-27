#include <iostream>				/* std::cout std::cerr */
#include <sstream>
#include <fstream>
#include "StorageBin.h"
#include "SS.h"
#ifndef boolean
typedef unsigned char boolean;
#endif
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Utils.h"


#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPTION_DEFAULT2 -5

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_line_LDBLEs(char *next_char, LDBLE ** d, int *count_d, int *count_alloc)
/* ---------------------------------------------------------------------- */
{
	int i, j, l, n;
	LDBLE value;
    std::string token;

	for (;;)
	{
        j = copy_token(token, next_char, &l);
		if (j == EMPTY)
			break;
		if (j != DIGIT)
			return (ERROR);
		if (replace("*", " ", token) == TRUE)
		{
//			if (sscanf(token, "%d" SCANFORMAT, &n, &value) != 2)
//				return (ERROR);
		}
		else
		{
//			sscanf(token, SCANFORMAT, &value);
			n = 1;
		}
		for (;;)
		{
			if ((*count_d) + n > (*count_alloc))
			{
				*count_alloc *= 2;
				*d = (LDBLE *)PHRQ_realloc(*d, (size_t)(*count_alloc) * sizeof(LDBLE));
				if (*d == NULL)
					malloc_error();
			}
			else
				break;
		}
		for (i = 0; i < n; i++)
			(*d)[(*count_d) + i] = value;
		*count_d += n;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
dump_cpp(void)
/* ---------------------------------------------------------------------- */
{
	/*
	* dumps solution compositions to file
	*/

	int l;

	if (dump_in == FALSE || pr.dump == FALSE)
		return (OK);

	cxxStorageBin phreeqcBin(phrq_io);
	phreeqc2cxxStorageBin(phreeqcBin);

	std::ofstream fs(dump_file_name_cpp.c_str());
	if (!fs.is_open())
	{
		error_string = sformatf("Can`t open file, %s.", dump_file_name_cpp.c_str());
		input_error++;
		error_msg(error_string, CONTINUE);
		return (OK);
	}

	fs << "# Dumpfile" << "\n" << "# Transport simulation " << simul_tr << "  Shift " << transport_step << "\n" << "#" << "\n";
	phreeqcBin.dump_raw(fs, 0);
	fs << "END" << "\n";

    std::string token;
//	sprintf(token, "KNOBS\n");
	fs << token;
//	sprintf(token, "\t-iter%15d\n", itmax);
	fs << token;
//	sprintf(token, "\t-tol %15.3e\n", (double)ineq_tol);
	fs << token;
//	sprintf(token, "\t-step%15.3e\n", (double)step_size);
	fs << token;
//	sprintf(token, "\t-pe_s%15.3e\n", (double)pe_step_size);
	fs << token;
//	sprintf(token, "\t-diag      ");
	fs << token;
	if (diagonal_scale == TRUE)
	{
//		sprintf(token, "true\n");
		fs << token;
	}
	else
	{
//		sprintf(token, "false\n");
		fs << token;
	}
	std::map < int, SelectedOutput >::iterator so_it = SelectedOutput_map.begin();
	for (; so_it != SelectedOutput_map.end(); so_it++)
	{
		current_selected_output = &(so_it->second);

//		sprintf(token, "SELECTED_OUTPUT %d\n", current_selected_output->Get_n_user());
		fs << token;
		//sprintf(token, "\t-file  %-15s\n", "sel_o$$$.prn");
		//fs << token;
		fs << "\t-file  " << "sel_o$$$" << current_selected_output->Get_n_user() << ".prn\n";
		//if (punch.count_totals != 0)
		if (current_selected_output->Get_totals().size() > 0)
		{
//			sprintf(token, "\t-tot ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_totals().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_totals()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_molalities().size() > 0)
		{
//			sprintf(token, "\t-mol ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_molalities().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_molalities()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_activities().size() > 0)
		{
//			sprintf(token, "\t-act ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_activities().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_activities()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_pure_phases().size() > 0)
		{
//			sprintf(token, "\t-equ ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_pure_phases().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_pure_phases()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_si().size() > 0)
		{
//			sprintf(token, "\t-si ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_si().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_si()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_gases().size() > 0)
		{
//			sprintf(token, "\t-gas ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_gases().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_gases()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_s_s().size() > 0)
		{
//			sprintf(token, "\t-solid_solutions ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_s_s().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_s_s()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_kinetics().size() > 0)
		{
//			sprintf(token, "\t-kin ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_kinetics().size(); i++)
			{
//				sprintf(token, "  %s", current_selected_output->Get_kinetics()[i].first.c_str());
				fs << token;
			}
//			sprintf(token, "\n");
			fs << token;
		}
	}
//	sprintf(token, "TRANSPORT\n");
	fs << token;
//	sprintf(token, "\t-cells %6d\n", count_cells);
	fs << token;
//	sprintf(token, "\t-shifts%6d%6d\n", count_shifts, ishift);
	fs << token;
//	sprintf(token, "\t-output_frequency %6d\n", print_modulus);
	fs << token;
//	sprintf(token, "\t-selected_output_frequency %6d\n",
//		punch_modulus);
	fs << token;
//	sprintf(token, "\t-bcon  %6d%6d\n", bcon_first, bcon_last);
	fs << token;
//	sprintf(token, "\t-timest %13.5e\n", (double)timest);
	fs << token;
	if (!high_precision)
	{
//		sprintf(token, "\t-diffc  %13.5e\n", (double)diffc);
		fs << token;
	}
	else
	{
//		sprintf(token, "\t-diffc  %20.12e\n", (double)diffc);
		fs << token;
	}
//	sprintf(token, "\t-tempr  %13.5e\n", (double)tempr);
	fs << token;
	if (correct_disp == TRUE)
	{
//		sprintf(token, "\t-correct_disp %s\n", "True");
		fs << token;
	}
	else
	{
//		sprintf(token, "\t-correct_disp %s\n", "False");
		fs << token;
	}
//	sprintf(token, "\t-length\n");
	fs << token;
	for (int i = 1; i <= count_cells; i++)
	{
//		sprintf(token, "%12.3e", (double)cell_data[i].length);
		fs << token;
		if (i > 0 && (i % 8) == 0)
		{
//			sprintf(token, "\n");
			fs << token;
		}
	}
//	sprintf(token, "\n");
	fs << token;
//	sprintf(token, "\t-disp\n");
	fs << token;
	for (int i = 1; i <= count_cells; i++)
	{
		if (!high_precision)
		{
//			sprintf(token, "%12.3e", (double)cell_data[i].disp);
			fs << token;
		}
		else
		{
//			sprintf(token, "%20.12e", (double)cell_data[i].disp);
			fs << token;
		}
		if (i > 0 && (i % 8) == 0)
		{
//			sprintf(token, "\n");
			fs << token;
		}
	}
//	sprintf(token, "\n");
	fs << token;
//	sprintf(token, "\t-punch_cells");
	fs << token;
	l = 0;
	for (int i = 0; i < all_cells; i++)
	{
		if (cell_data[i].punch != TRUE)
			continue;
//		sprintf(token, "  %d", i);
		fs << token;
		l++;
		if ((l % 20) == 0)
		{
//			sprintf(token, "\n");
			fs << token;
		}
	}
//	sprintf(token, "\n");
	fs << token;
//	sprintf(token, "\t-print_cells");
	fs << token;
	l = 0;
	for (int i = 0; i < all_cells; i++)
	{
		if (cell_data[i].print != TRUE)
			continue;
//		sprintf(token, "  %d", i);
		fs << token;
		l++;
		if ((l % 20) == 0)
		{
//			sprintf(token, "\n");
			fs << token;
		}
	}
//	sprintf(token, "\n");
	fs << token;
//	sprintf(token, "\t-dump            $$$.dmp\n");
	fs << token;
//	sprintf(token, "\t-dump_frequency  %d\n", dump_modulus);
	fs << token;
//	sprintf(token, "\t-dump_restart    %d\n", transport_step + 1);
	fs << token;

#if defined MULTICHART
	// user graphs
	chart_handler.dump(fs, 0);
#endif

//	sprintf(token, "END\n");
	fs << token;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
dump(void)
/* ---------------------------------------------------------------------- */
{
	/*
	* dumps solution compositions to file
	*/
	if (dump_in == FALSE || pr.dump == FALSE)
		return (OK);

	dump_cpp();
	return OK;

}
