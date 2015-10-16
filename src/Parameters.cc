#include "Parameters.h"

ParameterType::ParameterType() :
mShortSwitch(""), mLongSwitch(""), mShortDescription(""), mLongDescription(""), mTypeCode(STRING), mValue(""), mMinValue(""), mMaxValue(""), mNumStepsValue(""), mExponentialStepIncrease(false), mIsSet(false) {
}

void ParameterType::Parse(vector<string>& aParameterList) {
	string shortopt = "-" + mShortSwitch;
	string longopt = "--" + mLongSwitch;
	string shortopt_min = "-" + mShortSwitch + "-min";
	string shortopt_max = "-" + mShortSwitch + "-max";
	string longopt_min = "-" + mLongSwitch + "-min";
	string longopt_max = "-" + mLongSwitch + "-max";
	string shortopt_numsteps = "-" + mShortSwitch + "-num_steps";
	string longopt_numsteps = "-" + mLongSwitch + "-num_steps";

	string str_value = "";

	for (unsigned i = 0; i < aParameterList.size(); ++i) {
		if (aParameterList[i] == shortopt || aParameterList[i] == longopt || aParameterList[i] == shortopt_min || aParameterList[i] == shortopt_max || aParameterList[i] == longopt_min || aParameterList[i] == longopt_max || aParameterList[i] == shortopt_numsteps || aParameterList[i] == longopt_numsteps) {
			mIsSet = true;
			if (mTypeCode != FLAG) {
				if (i + 1 >= aParameterList.size())
					throw range_error("ERROR ParameterType::Parse: expected a value for option " + shortopt + "(" + longopt + ")");
				str_value = aParameterList[i + 1];
			}

			switch (mTypeCode) {
			case FLAG:
				str_value = "1";
				break;
			case REAL: {
				double value = strtod(str_value.c_str(), NULL);
				if (!value && str_value != "0")
					throw range_error("ERROR ParameterType::Parse: Value for option " + shortopt + " (" + longopt + ") must be of type real");
			}
			break;
			case INTEGER: {
				int value = strtol(str_value.c_str(), NULL, 10);
				if (!value && str_value != "0")
					throw range_error("ERROR ParameterType::Parse: Value for option " + shortopt + " (" + longopt + ") must be of type integer");
			}
			break;
			case POSITIVE_INTEGER: {
				unsigned value = strtoul(str_value.c_str(), NULL, 10);
				if (!value && str_value != "0")
					throw range_error("ERROR ParameterType::Parse: Value for option " + shortopt + " (" + longopt + ") must be of type positive integer");
			}
			break;
			case STRING: {
			}
			break;
			case LIST: {
				bool test = false;
				for (unsigned i = 0; i < mCloseValuesList.size(); ++i) {
					if (mCloseValuesList[i] == str_value) {
						test = true;
						break;
					}
				}
				if (test == false) {
					string value_list = "";
					for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
						value_list += mCloseValuesList[i] + " | ";
					value_list += mCloseValuesList.back();
					throw range_error("ERROR ParameterType::Parse: Value for option " + shortopt + " (" + longopt + ") must be one of: " + value_list + " instead of: " + mValue);
				}
			}
			break;
			default:
				throw range_error("ERROR ParameterType::Parse: Option " + shortopt + " (" + longopt + ") must be set to an appropriate value");
			}
		}
		if (aParameterList[i] == shortopt || aParameterList[i] == longopt) {
			mValue = str_value;
		}
		if (aParameterList[i] == shortopt_min || aParameterList[i] == longopt_min) {
			mMinValue = str_value;
		}
		if (aParameterList[i] == shortopt_max || aParameterList[i] == longopt_max) {
			mMaxValue = str_value;
		}
		if (aParameterList[i] == shortopt_numsteps || aParameterList[i] == longopt_numsteps) {
			mNumStepsValue = str_value;
		}
	}
}

void ParameterType::OutputCompact(ostream& out) const {
	if (mLongSwitch == "" && mShortSwitch == "")
		throw range_error("ERROR ParameterType::OutputCompact: Option does not have short nor long switch");
	if (mShortSwitch != "")
		out << "-" << mShortSwitch << " ";
	if (mLongSwitch != "")
		out << "--" << mLongSwitch << " ";
	if (mTypeCode == LIST) {
		string value_list = "";
		for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
			value_list += mCloseValuesList[i] + " | ";
		value_list += mCloseValuesList.back();
		out << value_list << " ";
	}
	if (mTypeCode != FLAG)
		if (mValue != "")
			out << "[" << mValue << "] ";
	if (mShortDescription != "")
		out << TAB << mShortDescription << " ";
	out << endl;
}

void ParameterType::OutputExtended(ostream& out) const {
	if (mLongSwitch == "" && mShortSwitch == "")
		throw range_error("ERROR ParameterType::OutputExtended: Option does not have short nor long switch");
	if (mShortSwitch != "")
		out << "-" << mShortSwitch << endl;
	if (mLongSwitch != "")
		out << "--" << mLongSwitch << endl;
	out << TAB << "Switch data type: ";
	switch (mTypeCode) {
	case FLAG:
		out << "flag" << endl;
		break;
	case REAL:
		out << "real number" << endl;
		break;
	case INTEGER:
		out << "integer number" << endl;
		break;
	case POSITIVE_INTEGER:
		out << "positive integer number" << endl;
		break;
	case STRING:
		out << "string" << endl;
		break;
	case LIST: {
		string value_list = "";
		for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
			value_list += mCloseValuesList[i] + " | ";
		value_list += mCloseValuesList.back();
		out << "One of the following " << mCloseValuesList.size() << " alternatives: " << value_list << endl;
	}
	break;
	default:
		throw range_error("ERROR ParameterType::OutputExtended: Unrecognized parameter type code [" + mTypeCode);
	}
	if (mTypeCode != FLAG) {
		if (mValue != "")
			out << TAB << "Default value: " << mValue << endl;
	}
	if (mLongDescription != "")
		out << TAB << "Description: " << mLongDescription << endl;
	else if (mShortDescription != "")
		out << TAB << "Description: " << mShortDescription << endl;
	out << endl;
}

//------------------------------------------------------------------------------------------------------
Parameters::Parameters() {
	SetupOptions();
}

void Parameters::SetupOptions() {
	{
		ParameterType param;
		param.mShortSwitch = "h";
		param.mLongSwitch = "help";
		param.mShortDescription = "Prints compact help.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
	}
	{
		ParameterType param;
		param.mShortSwitch = "H";
		param.mLongSwitch = "Help";
		param.mShortDescription = "Prints extended help.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
	}
	{
		ParameterType param;
		param.mShortSwitch = "a";
		param.mLongSwitch = "action";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "";
		param.mCloseValuesList.push_back("CLUSTER");
		param.mCloseValuesList.push_back("CLASSIFY");
		param.mCloseValuesList.push_back("TEST");


		mOptionList.insert(make_pair(param.mLongSwitch, param));

		mActionOptionList.insert(make_pair(CLUSTER, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(CLASSIFY, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TEST, vector<ParameterType*>()));

		string txt;
		txt = "Neighborhood Subgraph Pairwise Decomposition Kernel see: Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010.";
		mActionReferences.insert(make_pair(CLUSTER, txt));
		mActionReferences.insert(make_pair(CLASSIFY, txt));
		mActionReferences.insert(make_pair(TEST, txt));
		//Summaries
		txt = "Extract explicit feature representation using graph kernel decomposition.\n"
				"And nearest neighbors are efficiently identified with a locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(CLUSTER, txt));
		txt = "Clustering/classification of genomic and metagenomic sequences with locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(CLASSIFY, txt));
		txt = "Clustering/classification test for development.";
		mActionSummary.insert(make_pair(TEST, txt));
	}
	{
		ParameterType param;
		param.mShortSwitch = "y";
		param.mLongSwitch = "output_directory_path";
		param.mShortDescription = "";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "i";
		param.mLongSwitch = "input_data_file_name";
		param.mShortDescription = "";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "f";
		param.mLongSwitch = "file_type";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "FASTA";
		param.mCloseValuesList.push_back("STRINGSEQ");
		param.mCloseValuesList.push_back("FASTA");

		mOptionList.insert(make_pair(param.mLongSwitch, param));

		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}

	{
		ParameterType param;
		param.mShortSwitch = "k";
		param.mLongSwitch = "kernel_type";
		param.mLongDescription = "\n\n"
				"STRING=Neighborhood Subgraph Pairwise Decomposition Kernel for sequences - see: Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010.";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "STRING";
		param.mCloseValuesList.push_back("STRING");

		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "s";
		param.mLongSwitch = "suffix";
		param.mShortDescription = "Suffix string for all output files.";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "b";
		param.mLongSwitch = "hash_bit_size";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "30";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "r";
		param.mLongSwitch = "radius";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "2";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "d";
		param.mLongSwitch = "distance";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "5";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "V";
		param.mLongSwitch = "verbose";
		param.mShortDescription = "Outputs more information.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "cluster_type";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "APPROX_DENSE_CENTERS";
		param.mCloseValuesList.push_back("APPROX_DENSE_CENTERS");
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "F";
		param.mLongSwitch = "num_hash_functions";
		param.mShortDescription = "Number of random hash functions and length of signature";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "100";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "pure_approximate_sim";
		param.mShortDescription = "Calculation of the approximate neighborhood of an instance is only based on minHash signatures. Parameter value is between ]0..1] and determines the minimal fraction of common min-hash bins/hashes (similarity) of elements in the approximate neighborhood. 0=OFF and higher values also turn off e.g. -x etc.";
		param.mTypeCode = REAL;
		param.mValue = "0.5";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_repeat_hash_functions";
		param.mShortDescription = "Further split up feature space and repeat MinHash <VALUE> times for each slot. 1 is normal MinHash, 0 uses maximal repeated hash functions (equal to max=NumHashFunctions*NumHashShingles)";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "z";
		param.mLongSwitch = "max_size_bin";
		param.mShortDescription = "Maximum number of instances per bin in reverse index. When a bin contains references to more instances than this quantity, then the bin is erased. The ratio is that in these case the corresponding feature is common to too many instances and it is therefore not informative. Moreover the run-times become non sub-linear if each approximate neighborhood query implies to check a significant fraction of the data set.";
		param.mTypeCode = REAL;
		param.mValue = "1000";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_hash_shingles";
		param.mShortDescription = "Additional compression of an instance signature. Number of signature hash values combined into one hash values. Useful to identify highly similar instances.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	//--------------------------------------
	{
		ParameterType param;
		param.mShortSwitch = "R";
		param.mLongSwitch = "random_seed";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "numThreads";
		param.mShortDescription = "Used number of threads (OpenMP/std::thread) - 0 is max. hardware concurrency";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "seq_window";
		param.mShortDescription = "For FASTA files only! Defines the fragment length of sequences used for encoding! 0 = use full always the full length sequence.";
		param.mTypeCode = POSITIVE_INTEGER;
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			param.mValue = "100";
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "seq_clip";
		param.mShortDescription = "Clip this number of NT from each side of the seq, only applies if seq_window = 0, useful for seqs for classification to match seq_window of index";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "seq_shift";
		param.mShortDescription = "For FASTA files only! Defines the fraction as length of seq_window to generate sequence fragments used for classification! 0..1 Can be used to have a different shift for classification and indexing";
		param.mTypeCode = REAL;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "index_seq_shift";
		param.mShortDescription = "For FASTA files only! Defines the fraction as length of seq_window to generate sequence fragments used for encoding! 0..1";
		param.mTypeCode = REAL;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}

	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "min_radius";
		param.mShortDescription = "STRING kernel only: Lower bound radius used for feature encoding";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "min_distance";
		param.mShortDescription = "STRING kernel only: Lower bound distance used for feature encoding";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "index_bed";
		param.mShortDescription = "MinHash histogram reverse index is build from this data; The 4th col of BED entry denotes the bin in the histogram";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "index_seqs";
		param.mShortDescription = "MinHash histogram reverse index is build from this data; contains the corresponding seqs for regions given in BED file (--index_bed)";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "dense_center_names_file";
		param.mShortDescription = "list with seq names used as cluster centers; names have to match either fasta header (file_type FASTA) or idx (file_type STRINGSEQ); note: a dense cluster can contain not listed seqs";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "no_index_cache_file";
		param.mShortDescription = "Do not save/read MinHash reverse index to/from <index_data_file>.bhi file.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLASSIFY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "write_approx_neighbors";
		param.mShortDescription = "Writes all found approximate neighbors for each instance to file 'approx_dense_centers>'.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
}

void Parameters::Usage(string aCommandName, string aCompactOrExtended) {
	cout << SEP << endl << PROG_NAME << endl << "Version: " << PROG_VERSION << endl << "Last Update: " << PROG_DATE << endl << PROG_CREDIT << endl << SEP << endl;
	cout << "-v [--version] outputs current program version" << endl << endl;
	if (mActionSummary[mActionCode] != "")
		cout << "SUMMARY:" << endl << mActionSummary[mActionCode] << endl << SEP << endl;

	if (mOptionList["action"].mIsSet == false) {
		cout << "OPTIONS:" << endl;

		if (aCompactOrExtended == "EXTENDED")
			mOptionList["action"].OutputExtended(cout);
		else
			mOptionList["action"].OutputCompact(cout);
	} else {
		cout << "ACTION: " << mAction << endl;
		cout << "OPTIONS:" << endl;

		vector<ParameterType*> & vec = mActionOptionList[mActionCode];
		for (unsigned i = 0; i < vec.size(); ++i) {
			assert(vec[i]);
			ParameterType& param = (*vec[i]);
			if (aCompactOrExtended == "EXTENDED")
				param.OutputExtended(cout);
			else
				param.OutputCompact(cout);
		}
		cout << SEP << endl << "REFERENCES:" << endl << mActionReferences[mActionCode] << endl << SEP << endl;
	}
	exit(0);
}

void Parameters::Init(int argc, const char** argv) {
	if (argc == 1) {
		cout << "Use -h for compact help and -H for extended help." << endl;
		exit(1);
	}

	//convert argc in an option string vector
	vector<string> options;
	for (int i = 1; i < argc; i++)
		options.push_back(argv[i]);

	//parse the option string vector
	for (map<string, ParameterType>::iterator it = mOptionList.begin(); it != mOptionList.end(); ++it) {
		ParameterType& param = it->second;
		param.Parse(options);
	}

	//set the boolean parameters to a default value of false
	mVerbose = false;
	mNoIndexCacheFile = false;
	mWriteApproxNeighbors = false;
	//set the data members of Parameters according to user choice
	for (map<string, ParameterType>::iterator it = mOptionList.begin(); it != mOptionList.end(); ++it) {
		ParameterType& param = it->second;
		if (param.mIsSet) {
			if (param.mShortSwitch == "i")
				mInputDataFileName = param.mValue;
			if (param.mShortSwitch == "V")
				mVerbose = true;
			if (param.mLongSwitch == "write_approx_neighbors")
				mWriteApproxNeighbors = true;
			if (param.mLongSwitch == "no_index_cache_file")
				mNoIndexCacheFile = true;
		}


		if (param.mLongSwitch == "output_directory_path")
			mDirectoryPath = param.mValue;
		if (param.mShortSwitch == "a")
			mAction = param.mValue;
		if (param.mShortSwitch == "f")
			mFileType = param.mValue;
		if (param.mShortSwitch == "k")
			mKernelType = param.mValue;
		if (param.mShortSwitch == "s")
			mSuffix = param.mValue;
		if (param.mLongSwitch == "random_seed")
			mRandomSeed = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "r")
			mRadius = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "d")
			mDistance = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "b")
			mHashBitSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "F")
			mNumHashFunctions = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "z")
			mMaxSizeBin = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "fraction_center_scan")
			mFractionCenterScan = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "num_repeat_hash_functions")
			mNumRepeatsHashFunction = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "num_hash_shingles")
			mNumHashShingles = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "cluster_type")
			mClusterType = param.mValue;
		if (param.mLongSwitch == "numThreads")
			mNumThreads = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "max_fraction_of_dataset")
			mMaxFractionOfDataset = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "index_bed")
			mIndexBedFile = param.mValue;
		if (param.mLongSwitch == "index_seqs")
			mIndexSeqFile = param.mValue;
		if (param.mLongSwitch == "seq_shift")
			mSeqShift = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "index_seq_shift")
			mIndexSeqShift = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "seq_window")
			mSeqWindow = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "seq_clip")
			mSeqClip = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "min_radius")
			mMinRadius = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "min_distance")
			mMinDistance = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "pure_approximate_sim")
			mPureApproximateSim = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "dense_center_names_file")
			mDenseCenterNamesFile = param.mValue;
	}

	//convert action string to action code
	if (mAction == "")
		mActionCode = NULL_ACTION;
	else if (mAction == "CLUSTER")
		mActionCode = CLUSTER;
	else if (mAction == "CLASSIFY")
		mActionCode = CLASSIFY;
	else if (mAction == "TEST")
		mActionCode = TEST;
	else
		throw range_error("ERROR Parameters::Init: Unrecognized action: <" + mAction + ">");

	//convert file type string to file type code
	if (mFileType == "STRINGSEQ")
		mFileTypeCode = STRINGSEQ;
	else if (mFileType == "FASTA")
		mFileTypeCode = FASTA;
	else
		throw range_error("ERROR Parameters::Init: Unrecognized file type: <" + mFileType + ">");

	//check for help request
	for (unsigned i = 0; i < options.size(); ++i) {
		if (options[i] == "-h" || options[i] == "--help") {
			Usage(argv[0], "COMPACT");
			exit(1);
		}
		if (options[i] == "-H" || options[i] == "--Help") {
			Usage(argv[0], "EXTENDED");
			exit(1);
		}
		if (options[i] == "-v" || options[i] == "--version") {
			cout << "Version: " << PROG_VERSION << endl;
			exit(1);
		}
	}

	//check that set parameters are compatible
	if (mInputDataFileName == "")
		throw range_error("ERROR Parameters::Init: -i <input data file name> is missing.");
}
