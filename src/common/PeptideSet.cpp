#include "PeptideSet.h"

void PeptideSet::read(std::filesystem::path path)
{
	std::ifstream input(path.string());

	if (!input.is_open())
	{
		throw std::runtime_error("Could not open input");
	}

	std::string line, id, sequence;
	bool expect_sequence = false, expect_id = true;
	uint64_t line_count = 0;

	while (!input.eof())
	{
		std::getline(input, line); 
		if (line.empty() || isspace(line[0])) continue;
		++line_count;

		switch (line[0])
		{
		case ';': break;
		case '>':
			if (expect_sequence) throw std::runtime_error("expected sequence in the fasta file, got another ID");
			expect_sequence = true;
			expect_id = false;
			id = line.substr(1);
			break;
		default:
			expect_sequence = false;
			sequence = line;
			if (expect_id)
			{
				expect_id = false;
				id = std::string("L") + std::to_string(line_count);
			}

			emplace_back( sequence, id, (short)('f' - 'a') );
		}
	}
}

void PeptideSet::write(std::string_view path)
{
	std::ofstream output(path.data());
	for (auto&& peptide : *this)
	{
		output << ">" << peptide.id << std::endl;
		output << peptide.sequence << std::endl;
	}
}
