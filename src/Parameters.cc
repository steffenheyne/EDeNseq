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
		param.mCloseValuesList.push_back("TRAIN");
		param.mCloseValuesList.push_back("TEST");
		param.mCloseValuesList.push_back("TEST_PART");
		param.mCloseValuesList.push_back("CROSS_VALIDATION");
		param.mCloseValuesList.push_back("BIAS_VARIANCE_DECOMPOSITION");
		param.mCloseValuesList.push_back("PARAMETERS_OPTIMIZATION");
		param.mCloseValuesList.push_back("LEARNING_CURVE");
		param.mCloseValuesList.push_back("FEATURE");
		param.mCloseValuesList.push_back("FEATURE_PART");
		param.mCloseValuesList.push_back("FEATURE_SCALED");
		param.mCloseValuesList.push_back("MATRIX");
		param.mCloseValuesList.push_back("EMBED");
		param.mCloseValuesList.push_back("TARGET_ALIGNMENT");
		param.mCloseValuesList.push_back("CLUSTER");
		param.mCloseValuesList.push_back("NEAREST_NEIGHBOR");
		param.mCloseValuesList.push_back("MIN_HASH");
		param.mCloseValuesList.push_back("SEMI_SUPERVISED");
		param.mCloseValuesList.push_back("PRE_CLUSTER_FILTER");
		param.mCloseValuesList.push_back("DENSITY");
		param.mCloseValuesList.push_back("TRAIN_MULTILABEL");
		param.mCloseValuesList.push_back("TEST_MULTILABEL");
		param.mCloseValuesList.push_back("GENOME");


		mOptionList.insert(make_pair(param.mLongSwitch, param));

		mActionOptionList.insert(make_pair(TRAIN, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TEST, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TEST_PART, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(CROSS_VALIDATION, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(BIAS_VARIANCE_DECOMPOSITION, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(LEARNING_CURVE, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(PARAMETERS_OPTIMIZATION, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(FEATURE, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(FEATURE_PART, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(FEATURE_SCALED, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(MATRIX, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(EMBED, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TARGET_ALIGNMENT, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(CLUSTER, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(NEAREST_NEIGHBOR, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(MIN_HASH, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(SEMI_SUPERVISED, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(PRE_CLUSTER_FILTER, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(DENSITY, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TRAIN_MULTILABEL, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(TEST_MULTILABEL, vector<ParameterType*>()));
		mActionOptionList.insert(make_pair(GENOME, vector<ParameterType*>()));

		string txt;
		//References
		txt = "The code for Stochastic Gradient Descent SVM is adapted from http://leon.bottou.org/projects/sgd. Léon Bottou and Yann LeCun, ''Large Scale Online Learning'', Advances in Neural Information Processing Systems 16, Edited by Sebastian Thrun, Lawrence Saul and Bernhard Schölkopf, MIT Press, Cambridge, MA, 2004.\n";
		mActionReferences.insert(make_pair(TRAIN, txt));
		mActionReferences.insert(make_pair(TEST, txt));
		mActionReferences.insert(make_pair(TEST_PART, txt));
		mActionReferences.insert(make_pair(CROSS_VALIDATION, txt));
		txt = "The embedding method is adapted from L.Chen, A.Buja ''Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis'', Journal of the American Statistical Association, 2009.\n";
		mActionReferences.insert(make_pair(EMBED, txt));
		txt = "The semi supervised treatement is adapted from D. Zhou, J. Weston, A. Gretton, O. Bousquet and B. Schoelkopf, ''Ranking on Data Manifolds'', NIPS 2005 ";
		mActionReferences.insert(make_pair(SEMI_SUPERVISED, txt));

		//Summaries
		txt = "The linear model is induced using the accelerated stochastic gradient descent technique by Léon Bottou and Yann LeCun.\n"
				"When the target information is 0, a self-training algorithm is used to impute a positive or negative class to the unsupervised instances.\n"
				"If the target information is imbalanced a minority class oversampling technique is used to rebalance the training set.";
		mActionSummary.insert(make_pair(TRAIN, txt));
		mActionSummary.insert(make_pair(TEST, txt));
		mActionSummary.insert(make_pair(TEST_PART, txt));
		mActionSummary.insert(make_pair(CROSS_VALIDATION, txt));
		mActionSummary.insert(make_pair(BIAS_VARIANCE_DECOMPOSITION, txt));
		mActionSummary.insert(make_pair(LEARNING_CURVE, txt));
		txt = "Target information is iteratively spread to unsupervised nearest neighbors in a discounted fashion (using the kernel score).\n"
				"Nearest neighbors are efficiently identified with a locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(SEMI_SUPERVISED, txt));
		txt = "A multidimensional scaling technique is employed that treats nearest neighbors and distant instances differently. Distant instances are sampled randomly and are repelled with a uniform large force rather than according to the exact distance.\n"
				"Nearest neighbors are efficiently identified with a locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(EMBED, txt));
		txt = "Nearest neighbors are efficiently identified with a locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(NEAREST_NEIGHBOR, txt));
		txt = "Extract explicit feature representation using graph kernel decomposition.\n"
				"Using --smooth flag the features from the k-nearest neighbors are added, scaled by the corresponding kernel value.";
		mActionSummary.insert(make_pair(FEATURE, txt));
		txt = "Clustering/classification of genomic and metagenomic sequences with locality sensitive hashing technique.";
		mActionSummary.insert(make_pair(GENOME, txt));
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "t";
		param.mLongSwitch = "target_file_name";
		param.mShortDescription = "";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "extended_matrix_information";
		param.mShortDescription = "In MATRIX.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "A";
		param.mLongSwitch = "row_index_file_name";
		param.mShortDescription = "In MATRIX and in NEAREST_NEIGHBOR.";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "B";
		param.mLongSwitch = "col_index_file_name";
		param.mShortDescription = "In MATRIX and in NEAREST_NEIGHBOR.";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "smooth";
		param.mShortDescription = "Adds rescaled features from nearest neighbors.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}

	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "smoother_param";
		param.mShortDescription = "Features from neighbors are scaled by the kernel value to the power value assigned to this switch.";
		param.mTypeCode = REAL;
		param.mValue = ".925";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
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
		param.mValue = "GRAPH";
		param.mCloseValuesList.push_back("SPARSE_VECTOR");
		param.mCloseValuesList.push_back("GRAPH");
#ifdef USEOBABEL
		param.mCloseValuesList.push_back("MOLECULAR_GRAPH");
#endif
		param.mCloseValuesList.push_back("SEQUENCE");
		param.mCloseValuesList.push_back("STRINGSEQ");
		param.mCloseValuesList.push_back("FASTA");

		mOptionList.insert(make_pair(param.mLongSwitch, param));

		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}

	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "binary_file_type";
		param.mShortDescription = "";
		param.mTypeCode = FLAG;
		param.mValue = "0";

		mOptionList.insert(make_pair(param.mLongSwitch, param));

		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "k";
		param.mLongSwitch = "kernel_type";
		param.mLongDescription = "\n\n"
				"NSPDK=Neighborhood Subgraph Pairwise Decomposition Kernel see: Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010."
				"\n\n"
				"WDK=Weighted Decomposition Kernel, see: Sauro Menchetti, Fabrizio Costa, Paolo Frasconi, '' ''"
				"\n\n"
				"PBK=Poly Bias Kernel"
				"\n\n"
				"USPK="
				"\n\n"
				"DDK=Decompositional DAG graph Kernel, see: G. Da San Martino, N. Navarin and A. Sperduti , ''A tree-based kernel for graphs'', Proceedings of the Twelfth SIAM International Conference on Data Mining, Anaheim, California, April 26 - 28, 2012, p. 975-986."
				"\n\n"
				"NSDDK="
				"\n\n"
				"ANSDDK="
				"\n\n"
				"SK=";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "NSPDK";
		param.mCloseValuesList.push_back("NSPDK");
		param.mCloseValuesList.push_back("WDK");
		param.mCloseValuesList.push_back("PBK");
		param.mCloseValuesList.push_back("USPK");
		param.mCloseValuesList.push_back("DDK");
		param.mCloseValuesList.push_back("NSDDK");
		param.mCloseValuesList.push_back("ANSDDK");
		param.mCloseValuesList.push_back("SK");
		param.mCloseValuesList.push_back("STRING");

		mOptionList.insert(make_pair(param.mLongSwitch, param));

		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "g";
		param.mLongSwitch = "graph_type";
		param.mShortDescription = "";
		param.mTypeCode = LIST;
		param.mValue = "UNDIRECTED";
		param.mCloseValuesList.push_back("DIRECTED");
		param.mCloseValuesList.push_back("UNDIRECTED");
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "use_real_vector_information";
		param.mShortDescription = "";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_random_projections";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "no_normalization";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "min_kernel";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
#ifdef USEOBABEL
	{
		ParameterType param;
		param.mShortSwitch = "o";
		param.mLongSwitch = "open_babel_file_format";
		param.mShortDescription = "MOLECULAR_GRAPH parameter.";
		param.mTypeCode = STRING;
		param.mValue = "sdf";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
#endif
	{
		ParameterType param;
		param.mShortSwitch = "m";
		param.mLongSwitch = "model_file_name";
		param.mShortDescription = "";
		param.mTypeCode = STRING;
		param.mValue = "model";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
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
		param.mValue = "15";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "vertex_degree_threshold";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "7";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "M";
		param.mLongSwitch = "sequence_degree";
		param.mShortDescription = "SEQUENCE data type.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "sequence_token";
		param.mShortDescription = "Labels are strings separated by spaces.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "sequence_multi_line";
		param.mShortDescription = "The annotation is encoded on subsequent rows.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "sequence_pairwise_interaction";
		param.mShortDescription = "Abstraction vertices relating all pairs of disjointed sequences are inserted.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "l";
		param.mLongSwitch = "lambda";
		param.mShortDescription = "Stochastic gradient descend algorithm.";
		param.mTypeCode = REAL;
		param.mValue = "1e-4";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "e";
		param.mLongSwitch = "epochs";
		param.mShortDescription = "Stochastic gradient descend algorithm.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN_MULTILABEL];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "O";
		param.mLongSwitch = "sparsification_num_iterations";
		param.mShortDescription = "In training.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "C";
		param.mLongSwitch = "topological_regularization_num_neighbors";
		param.mShortDescription = "In training.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "L";
		param.mLongSwitch = "topological_regularization_decay_rate";
		param.mShortDescription = "In training.";
		param.mTypeCode = REAL;
		param.mValue = "0.01";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "c";
		param.mLongSwitch = "num_cross_validation_folds";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "p";
		param.mLongSwitch = "num_evaluation_points";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "Y";
		param.mLongSwitch = "num_line_search_iterations";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "3";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
	}
	{
		ParameterType param;
		param.mShortSwitch = "S";
		param.mLongSwitch = "num_iterations";
		param.mShortDescription = "In semi-supervised setting.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "3";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "semi_supervised_alpha";
		param.mShortDescription = "In SEMI_SUPERVISED";
		param.mTypeCode = REAL;
		param.mValue = ".99";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "T";
		param.mLongSwitch = "threshold";
		param.mShortDescription = "In semi-supervised setting. Only the top and low quantile will be used as positives and negative instances. A threshold of 1 means that all unsupervised instaces are used in the next phase.";
		param.mTypeCode = REAL;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "only_positive";
		param.mShortDescription = "In semi-supervised setting. Induce only positive class instances.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "only_negative";
		param.mShortDescription = "In semi-supervised setting. Induce only negative class instances.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "N";
		param.mLongSwitch = "num_of_random_restarts";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "I";
		param.mLongSwitch = "num_of_iterations";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "5000";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "D";
		param.mLongSwitch = "dimensionality";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "2";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "E";
		param.mLongSwitch = "epsilon";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = REAL;
		param.mValue = "0.01";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "K";
		param.mLongSwitch = "neighborhood_size";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "refine_neighborhood_size";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "n";
		param.mLongSwitch = "non_neighborhood_size";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "100";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "G";
		param.mLongSwitch = "neighborhood_size_range";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "U";
		param.mLongSwitch = "tau";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = REAL;
		param.mValue = "0.0005";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "u";
		param.mLongSwitch = "tau_exponent_range";
		param.mShortDescription = "In EMBED.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "1";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "V";
		param.mLongSwitch = "verbose";
		param.mShortDescription = "Outputs the graphs and the feature encodings.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "W";
		param.mLongSwitch = "minimal_output";
		param.mShortDescription = "Output only results.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
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
		param.mValue = "K_QUICK_SHIFT";
		param.mCloseValuesList.push_back("K_QUICK_SHIFT");
		param.mCloseValuesList.push_back("DENSE_CENTERS");
		param.mCloseValuesList.push_back("DENSE_CONNECTED_CENTERS");
		param.mCloseValuesList.push_back("APPROXIMATION_ACCURACY");
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
		param.mShortSwitch = "";
		param.mLongSwitch = "approximate";
		param.mShortDescription = "Use approximate nearest neighbor queries. Using hashing techniques the complexity of these queries is constant time, rather than linear.";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "F";
		param.mLongSwitch = "num_hash_functions";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "400";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_repeat_hash_functions";
		param.mShortDescription = "In approximate neighborhood";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "z";
		param.mLongSwitch = "max_size_bin";
		param.mShortDescription = "In CLUSTER.  Expressed as the maximum number of instances per bin. When a bin contains references to more instances than this quantity, then the bin is erased. The ratio is that in these case the corresponding feature is common to too many instances and it is therefore not informative. Moreover the run-times become non sub-linear if each approximate neighborhood query implies to check a significant fraction of the data set.";
		param.mTypeCode = REAL;
		param.mValue = "1000";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "eccess_neighbor_size_factor";
		param.mShortDescription = "In CLUSTER. Expressed as a multiplicative factor w.r.t. the neighborhood size required. It means that the approximate neighborhood query stops at the X most frequent instances, where X= eccess_neighbor_size_factor * neighborhood size.";
		param.mTypeCode = REAL;
		param.mValue = "5.0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "Q";
		param.mLongSwitch = "sample_size";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "20";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "x";
		param.mLongSwitch = "num_nearest_neighbors";
		param.mShortDescription = "In CLUSTER and in NEAREST_NEIGHBOR";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "X";
		param.mLongSwitch = "shared_neighborhood";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "fraction_center_scan";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = REAL;
		param.mValue = "0.5";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "w";
		param.mLongSwitch = "neighborhood_intersection_size";
		param.mShortDescription = "";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "3";
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
		param.mLongSwitch = "cluster_threshold";
		param.mShortDescription = "Maximum factor for parent similarity w.r.t. average neighbors similarity (normalized). A smaller value (e.g. 0.5) implies a stricter similarity criterion. In CLUSTER";
		param.mTypeCode = REAL;
		param.mValue = "1";
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
		param.mLongSwitch = "no_neighborhood_cache";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "no_minhash_cache";
		param.mShortDescription = "In CLUSTER";
		param.mTypeCode = FLAG;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[SEMI_SUPERVISED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "pure_approximate_sim";
		param.mShortDescription = "In CLUSTER. Calculation of the approximate neighborhood of an instance is only based on minHash signatures. Parameter value is between ]0..1] and determines the minimal fraction of common min-hash bins/hashes (similarity) of elements in the approximate neighborhood. 0=OFF and higher values also turn off e.g. -x etc.";
		param.mTypeCode = REAL;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_hash_shinglets";
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
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	//DD kernel family additional parameters
	{
		ParameterType param;
		param.mLongSwitch = "tree_lambda";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = REAL;
		param.mValue = "1.2";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mLongSwitch = "radius_two";
		param.mShortDescription = "Kernel parameter.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "2";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TEST_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_PART];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[FEATURE_SCALED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[MATRIX];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[TARGET_ALIGNMENT];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
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
			vector<ParameterType*>& vec = mActionOptionList[TRAIN];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[CROSS_VALIDATION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[BIAS_VARIANCE_DECOMPOSITION];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[LEARNING_CURVE];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "initial_embedding_coordinates_file_name";
		param.mShortDescription = "";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[EMBED];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_output_instances";
		param.mShortDescription = "Number of feature vectors to return after filtering using an efficient proxy to assess the probability for each instance to be part of a dense cluster.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "10000";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "size_random_dataset";
		param.mShortDescription = "Number of instances randomly generated to estimate the empirical probability density function.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "5000";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[DENSITY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "num_reshuffles";
		param.mShortDescription = "Number of times the density estimate procedure is performed.";
		param.mTypeCode = POSITIVE_INTEGER;
		param.mValue = "30";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[DENSITY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "max_fraction_of_dataset";
		param.mShortDescription = "Max fraction of dataset size used in density estimate. Its has a smoothing effect on the estimate.";
		param.mTypeCode = REAL;
		param.mValue = ".9";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[DENSITY];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "std_dev_random_sampling";
		param.mShortDescription = "Standard deviation for generator of random instances";
		param.mTypeCode = REAL;
		param.mValue = ".001";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[DENSITY];
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
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
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
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "seq_shift";
		param.mShortDescription = "For FASTA files only! Defines the fraction as length of seq_window to generate sequence fragments used for encoding! 0..1";
		param.mTypeCode = REAL;
		param.mValue = "0";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
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
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
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
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
	}
	{
		ParameterType param;
		param.mShortSwitch = "";
		param.mLongSwitch = "index_data_file";
		param.mShortDescription = "MinHash histogram reverse index is build from this data; <IDX> <FILENAME> per line";
		param.mTypeCode = STRING;
		param.mValue = "";
		mOptionList.insert(make_pair(param.mLongSwitch, param));
		{
			vector<ParameterType*>& vec = mActionOptionList[CLUSTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[NEAREST_NEIGHBOR];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[PRE_CLUSTER_FILTER];
			ParameterType& p = mOptionList[param.mLongSwitch];
			vec.push_back(&p);
		}
		{
			vector<ParameterType*>& vec = mActionOptionList[GENOME];
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
	mKernelNoNormalization = false;
	mMinKernel = false;
	mSemiSupervisedInduceOnlyPositive = false;
	mSemiSupervisedInduceOnlyNegative = false;
	mVerbose = false;
	mSequenceToken = false;
	mSequenceMultiLine = false;
	mSequencePairwiseInteraction = false;
	mMinimalOutput = false;
	mSharedNeighborhood = false;
	mBinaryFormat = false;
	mUseApproximate = false;
	mNoNeighborhoodCache = false;
	mNoMinHashCache = false;
	mSmooth = false;
	mExtendedMatrixInformation = false;
	mUseRealVectorInformation = false;

	//set the data members of Parameters according to user choice
	for (map<string, ParameterType>::iterator it = mOptionList.begin(); it != mOptionList.end(); ++it) {
		ParameterType& param = it->second;
		if (param.mIsSet) {
			if (param.mShortSwitch == "i")
				mInputDataFileName = param.mValue;
			if (param.mShortSwitch == "t")
				mTargetFileName = param.mValue;
			if (param.mLongSwitch == "no_normalization")
				mKernelNoNormalization = true;
			if (param.mLongSwitch == "min_kernel")
				mMinKernel = true;
			if (param.mLongSwitch == "only_positive")
				mSemiSupervisedInduceOnlyPositive = true;
			if (param.mLongSwitch == "only_negative")
				mSemiSupervisedInduceOnlyNegative = true;
			if (param.mShortSwitch == "V")
				mVerbose = true;
			if (param.mShortSwitch == "W")
				mMinimalOutput = true;
			if (param.mLongSwitch == "sequence_token")
				mSequenceToken = true;
			if (param.mLongSwitch == "sequence_multi_line")
				mSequenceMultiLine = true;
			if (param.mLongSwitch == "sequence_pairwise_interaction")
				mSequencePairwiseInteraction = true;
			if (param.mShortSwitch == "X")
				mSharedNeighborhood = true;
			if (param.mLongSwitch == "binary_file_type")
				mBinaryFormat = true;
			if (param.mLongSwitch == "approximate")
				mUseApproximate = true;
			if (param.mLongSwitch == "no_neighborhood_cache")
				mNoNeighborhoodCache = true;
			if (param.mLongSwitch == "no_minhash_cache")
				mNoMinHashCache = true;
			if (param.mLongSwitch == "smooth")
				mSmooth = true;
			if (param.mLongSwitch == "extended_matrix_information")
				mExtendedMatrixInformation = true;
			if (param.mLongSwitch == "use_real_vector_information")
				mUseRealVectorInformation = true;
		}


		if (param.mLongSwitch == "initial_embedding_coordinates_file_name")
			mEmbedFileName = param.mValue;
		if (param.mLongSwitch == "output_directory_path")
			mDirectoryPath = param.mValue;
		if (param.mShortSwitch == "a")
			mAction = param.mValue;
		if (param.mShortSwitch == "A")
			mRowIndexFileName = param.mValue;
		if (param.mShortSwitch == "B")
			mColIndexFileName = param.mValue;
		if (param.mShortSwitch == "m")
			mModelFileName = param.mValue;
		if (param.mShortSwitch == "f")
			mFileType = param.mValue;
		if (param.mShortSwitch == "k")
			mKernelType = param.mValue;
		if (param.mShortSwitch == "g")
			mGraphType = param.mValue;
		if (param.mShortSwitch == "o")
			mOpenBabelFormat = param.mValue;
		if (param.mShortSwitch == "s")
			mSuffix = param.mValue;
		if (param.mShortSwitch == "r")
			mRadius = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "d")
			mDistance = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "vertex_degree_threshold")
			mVertexDegreeThreshold = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "M")
			mSequenceDegree = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "b")
			mHashBitSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "l")
			mLambda = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "e")
			mEpochs = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "c")
			mCrossValidationNumFolds = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "p")
			mNumPoints = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "R")
			mRandomSeed = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "S")
			mSemiSupervisedNumIterations = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "T")
			mSemiSupervisedThreshold = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "N")
			mLMDSNumRandomRestarts = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "I")
			mLMDSNumIterations = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "D")
			mLMDSDimensionality = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "E")
			mLMDSIterationEpsilon = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "K")
			mLMDSNeighborhoodSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "n")
			mLMDSNonNeighborhoodSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "G")
			mLMDSNeighborhoodSizeRange = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "U")
			mLMDSTau = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "u")
			mLMDSTauExponentRange = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "O")
			mSparsificationNumIterations = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "C")
			mTopologicalRegularizationNumNeighbors = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "L")
			mTopologicalRegularizationRate = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "F")
			mNumHashFunctions = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "z")
			mMaxSizeBin = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "eccess_neighbor_size_factor")
			mEccessNeighborSizeFactor = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "Q")
			mSampleSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "x")
			mNumNearestNeighbors = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "fraction_center_scan")
			mFractionCenterScan = stream_cast<double>(param.mValue);
		if (param.mShortSwitch == "w")
			mNeighborhoodIntersectionSize = stream_cast<unsigned>(param.mValue);
		if (param.mShortSwitch == "Y")
			mNumLineSearchIterations = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "num_repeat_hash_functions")
			mNumRepeatsHashFunction = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "cluster_threshold")
			mClusterThreshold = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "pure_approximate_sim")
			mPureApproximateSim = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "num_hash_shinglets")
			mNumHashShinglets = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "cluster_type")
			mClusterType = param.mValue;
		if (param.mLongSwitch == "numThreads")
			mNumThreads = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "semi_supervised_alpha")
			mSemiSupervisedAlpha = stream_cast<double>(param.mValue);
		//DD kernel family additional parameters
		if (param.mLongSwitch == "tree_lambda")
			mTreeLambda = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "radius_two")
			mRadiusTwo = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "smoother_param")
			mSmootherParam = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "refine_neighborhood_size")
			mLMDSRefineNumNearestNeighbor = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "num_output_instances")
			mNumOutputInstances = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "size_random_dataset")
			mSizeRandomDataset = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "num_reshuffles")
			mNumReshuffles = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "max_fraction_of_dataset")
			mMaxFractionOfDataset = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "std_dev_random_sampling")
			mStdDevRandomSampling = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "num_random_projections")
			mNumRandomProjections = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "index_data_file")
			mIndexDataList = param.mValue;
		if (param.mLongSwitch == "seq_shift")
			mSeqShift = stream_cast<double>(param.mValue);
		if (param.mLongSwitch == "seq_window")
			mSeqWindow = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "min_radius")
			mMinRadius = stream_cast<unsigned>(param.mValue);
		if (param.mLongSwitch == "min_distance")
			mMinDistance = stream_cast<unsigned>(param.mValue);

	}

	//convert action string to action code
	if (mAction == "")
		mActionCode = NULL_ACTION;
	else if (mAction == "TRAIN")
		mActionCode = TRAIN;
	else if (mAction == "TEST")
		mActionCode = TEST;
	else if (mAction == "CROSS_VALIDATION")
		mActionCode = CROSS_VALIDATION;
	else if (mAction == "BIAS_VARIANCE_DECOMPOSITION")
		mActionCode = BIAS_VARIANCE_DECOMPOSITION;
	else if (mAction == "PARAMETERS_OPTIMIZATION")
		mActionCode = PARAMETERS_OPTIMIZATION;
	else if (mAction == "LEARNING_CURVE")
		mActionCode = LEARNING_CURVE;
	else if (mAction == "TEST_PART")
		mActionCode = TEST_PART;
	else if (mAction == "FEATURE")
		mActionCode = FEATURE;
	else if (mAction == "FEATURE_PART")
		mActionCode = FEATURE_PART;
	else if (mAction == "FEATURE_SCALED")
		mActionCode = FEATURE_SCALED;
	else if (mAction == "MATRIX")
		mActionCode = MATRIX;
	else if (mAction == "EMBED")
		mActionCode = EMBED;
	else if (mAction == "TARGET_ALIGNMENT")
		mActionCode = TARGET_ALIGNMENT;
	else if (mAction == "CLUSTER")
		mActionCode = CLUSTER;
	else if (mAction == "NEAREST_NEIGHBOR")
		mActionCode = NEAREST_NEIGHBOR;
	else if (mAction == "MIN_HASH")
		mActionCode = MIN_HASH;
	else if (mAction == "SEMI_SUPERVISED")
		mActionCode = SEMI_SUPERVISED;
	else if (mAction == "PRE_CLUSTER_FILTER")
		mActionCode = PRE_CLUSTER_FILTER;
	else if (mAction == "DENSITY")
		mActionCode = DENSITY;
	else if (mAction == "TRAIN_MULTILABEL")
		mActionCode = TRAIN_MULTILABEL;
	else if (mAction == "TEST_MULTILABEL")
		mActionCode = TEST_MULTILABEL;
	else if (mAction == "GENOME")
		mActionCode = GENOME;
	else
		throw range_error("ERROR Parameters::Init: Unrecognized action: <" + mAction + ">");

	//convert file type string to file type code
	if (mFileType == "GRAPH")
		mFileTypeCode = GRAPH;
	else if (mFileType == "SPARSE_VECTOR")
		mFileTypeCode = SPARSE_VECTOR;
#ifdef USEOBABEL
	else if (mFileType == "MOLECULAR_GRAPH")
		mFileTypeCode = MOLECULAR_GRAPH;
#endif
	else if (mFileType == "SEQUENCE")
		mFileTypeCode = SEQUENCE;
	else if (mFileType == "STRINGSEQ")
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

	if ((mAction == "TRAIN" || mAction == "CROSS_VALIDATION" || mAction == "LEARNING_CURVE" || mAction == "TARGET_ALIGNMENT" || mAction == "TRAIN_MULTILABEL") && mTargetFileName == "")
		throw range_error("ERROR Parameters::Init: -t <target file name> is missing.");

	if (mAction == "FEATURE" && mSmooth == true && (mRowIndexFileName == "" || mColIndexFileName == ""))
		throw range_error("ERROR Parameters::Init: When using --smooth the switches -A and -B  must be set.");

	if (mKernelType == "STRING" && mFileType != "STRINGSEQ" && mActionCode != CLUSTER)
		throw range_error("ERROR Parameters::Init: when using a STRING kernel the format must be STRINGSEQ");

	if (mKernelType != "STRING" && mFileType == "STRINGSEQ")
		throw range_error("ERROR Parameters::Init: when using the STRINGSEQ format the kernel must be of type STRING");

#ifndef USEOBABEL
	if (mFileType == "MOLECULAR_GRAPH")
		throw range_error("ERROR Parameters::Init: -f MOLECULAR_GRAPH is not enabled if compilation flag USEOBABEL is not set. Please recompile using -DUSEOBABEL");
#endif
}
