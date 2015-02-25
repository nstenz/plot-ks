#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use IO::Pipe;
use IO::Select;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Path qw(remove_tree);
use Storable qw(fd_retrieve store_fd);

# TODO: add support for multi-threading
my $previous_max_forks;
my $max_forks = get_free_cpus();

# Hash used to perform reverse translations
my %rev_codon_table = (
	S => qr/((AG[C|T])|(TC.))/,
	F => qr/TT[T|C]/,
	L => qr/((TT[A|G])|(CT.))/,
	Y => qr/TA[T|C]/,
	C => qr/TG[T|C]/,
	W => qr/TGG/,
	P => qr/CC./,
	H => qr/CA[T|C]/,
	Q => qr/CA[A|G]/,
	R => qr/((AG[A|G])|(CG.))/,
	I => qr/AT[T|C|A]/,
	M => qr/ATG/,
	T => qr/AC./,
	N => qr/AA[T|C]/,
	K => qr/AA[A|G]/,
	V => qr/GT./,
	A => qr/GC./,
	D => qr/GA[T|C]/,
	E => qr/GA[A|G]/,
	G => qr/GG./,
	X => qr/.../,
);

# Will contain PID of process running pipeline
my $pid;

# Allowed models 
my %models = (NG => 1, LWL => 1, LPB => 1, MLWL => 1, MLPB => 1, 
              GY => 1, YN => 1, MYN => 1, MS => 1, MA => 1);

# Default settings
my $model = "YN";
my $bin_size = 0.05;
my $match_length_threshold = 300; # Nucleotides

# Range for Ks plot
my $ks_min = 0;
my $ks_max = 3;
my $exclude_zero;

# Allow for reusing info from an old run
my $input_is_dir = 0;

# Name of output directory
my $project_name = "plot-ks-".time();

my $kaks_calc_in_dir = "kaks-in";
my $kaks_calc_out_dir = "kaks-out";
my $split_blat_out_dir = "blat-out";

# Check that dependencies can be found in user's PATH
my $r = check_path_for_exec("R");
my $blat = check_path_for_exec("blat");
my $transdecoder = check_path_for_exec("TransDecoder");
my $kaks_calculator = check_path_for_exec("KaKs_Calculator");

# Save how script was invoked
my $settings = "@ARGV";

# Parse command line options
GetOptions(
	"model|m:s" => \$model,
	"ks_min:f" => \$ks_min,
	"ks_max:f" => \$ks_max,
	"bin_size|b:f" => \$bin_size,
	"exclude_zero|x" => \$exclude_zero,
	"match_length_threshold|t:i" => \$match_length_threshold,
	"help|h" => \&help,
	"usage" => \&usage
);

# Error check input
my $transcriptome = shift;
die "You must specify a FASTA file containing transcriptome data for input OR the output directory from the previous execution of this script.\n".&usage if (!defined($transcriptome));
die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
die "Your Ks plot bin size can't be a negative number.\n".&usage if ($bin_size < 0);
die "Your Ks plot bin size is larger than your Ks range.\n".&usage if ($bin_size > $ks_max - $ks_min);
die "The model '$model' does not exist or is not usable by this script.\n".&usage if (!exists($models{$model}));

# Input is a previous run directory, reuse information
$input_is_dir++ if (-d $transcriptome);

my $input_root;
if (!$input_is_dir) {
	my $transcriptome_abs_path = abs_path($transcriptome);

	# Initialize our working directories and create a symlink to input
	mkdir($project_name) if (!-e $project_name) || die "Could not create directory '$project_name': $!.\n";

	mkdir("$project_name/$kaks_calc_in_dir") || die "Could not create directory '$kaks_calc_in_dir': $!.\n";
	mkdir("$project_name/$kaks_calc_out_dir") || die "Could not create directory '$kaks_calc_out_dir': $!.\n";
	mkdir("$project_name/$split_blat_out_dir") || die "Could not create directory '$split_blat_out_dir': $!.\n";

	$transcriptome =~ s/.*\///;
	system("ln -s $transcriptome_abs_path $ENV{PWD}/$project_name/$transcriptome");
	($input_root = $transcriptome) =~ s/(.*)\.\S+$/$1/;
}
else {
	$project_name = $transcriptome;
	my @contents = glob("$project_name/*");

	# Determine which transcriptome name by looking for a symlink
	my $found_name = 0;
	foreach my $file (@contents) {
		if (-l $file) {
			$file =~ s/\Q$project_name\E\///;
			$transcriptome = $file;
			$found_name = 1;
		}
	}
	($input_root = $transcriptome) =~ s/(.*)\.\S+$/$1/;

	# For simplicity we will reuse the same number of forks 
	if (-e "$project_name/$kaks_calc_in_dir") {
		@contents = glob("$project_name/*$kaks_calc_in_dir");
		$previous_max_forks = scalar(@contents);
	}

	die "Could not locate transcriptome in '$project_name'.\n" if (!$found_name);
}

chdir($project_name);

# Change how we handle Ctrl+C
$SIG{INT} = 'INT_handler';

# Set up output file names
my $blat_out = "self-blat.pslx";
#my $paralogs_seqs = "paralogs-t$match_length_threshold.atx";
#my $paralogs_ka_ks_out = "paralogs-t$match_length_threshold-m$model.kaks";
my $paralogs_seqs = "$kaks_calc_in_dir/paralogs-t$match_length_threshold.atx";
my $paralogs_ka_ks_out = "$kaks_calc_out_dir/paralogs-t$match_length_threshold-m$model.kaks";
my $final_ks_values = "paralogs-t$match_length_threshold-m$model.csv";

# Output plot name varies based on bin size and Ks range
my $ks_plot_name = "$input_root-t$match_length_threshold-m$model";
if ($exclude_zero) {
	if ($bin_size == 0) {
		$ks_plot_name .= "-density-range\(0-$ks_max].pdf";
	}
	else {
		$ks_plot_name .= "-b$bin_size-range\(0-$ks_max].pdf";
	}
}
else {
	if ($bin_size == 0) {
		$ks_plot_name .= "-density-range[0-$ks_max].pdf";
	}
	else {
		$ks_plot_name .= "-b$bin_size-range[0-$ks_max].pdf";
	}
}

$pid = fork();

# Child fork
if ($pid == 0) {

	# Cleans up if cancelled
	$SIG{TERM} = 'TERM_handler';

	# Echo script invocation
	logger("Invocation: perl plot-ks.pl $settings");

	# Translate transcriptome into most likely proteome
	run_transdecoder();

	# Identify paralogs
	run_self_blat();
	parse_self_blat_output();

	# Calculate Ks for paralogs
	run_kaks_calc();
	parse_ks_values();

	# Create the final plot
	create_ks_plot();
	exit(0);
}
else {
	waitpid($pid, 0);
}

sub run_transdecoder {

	my $transdecoder_out_dir = "transdecoder";

	# Check if we've already run transdecoder
	if (!-e "$transcriptome.transdecoder.pep" || !-e "$transcriptome.transdecoder.mRNA") {
		logger("\nRunning TransDecoder on '$transcriptome'...");
		system("$transdecoder -t $transcriptome --workdir $transdecoder_out_dir") && die;

		# Clean up unneeded files
		remove_tree($transdecoder_out_dir);

		logger("Completed TransDecoder.\n");
	}
	else {
		logger("\nUsing TransDecoder output from previous run.");
	}
}

sub run_self_blat {

	# Check if we've already run blat
	if (!-e $blat_out) {
		logger("Running self blat...");
		system("$blat $transcriptome.transdecoder.pep $transcriptome.transdecoder.pep -prot -out=pslx self-blat.pslx -noHead") && die;
		logger("Completed self blat.\n");
	}
	else {
		logger("Using blat output from previous run in '$blat_out'.");
	}
}

sub parse_self_blat_output_child {
	#my ($file, $align) = @_;
	my ($file, $handle, $align) = @_;

	my %align = %{$align};

	my $count = 0;

	#my %queries;
	my %matches;
	open(my $blat_out, "<", $file);
	while (my $line = <$blat_out>) {
		chomp($line);

		# Split the line on tabs and extract the information we want
		my @line = split("\t", $line);
		my ($query_name, $match_name, $query_align, $match_align) = 
			($line[9], $line[13], $line[21], $line[22]);

		# Don't want self matches
		next if ($query_name eq $match_name);
		
		# Reverse translate amino acids to their corresponding nucleotides

		my @query_align =  split(",", $query_align);
		my @match_align =  split(",", $match_align);

		my $query_nuc_align = $align{$query_name};
		my $match_nuc_align = $align{$match_name};

		my $trans_query_align;
		foreach my $align (@query_align) {
			$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
		}

		my $trans_match_align;
		foreach my $align (@match_align) {
			$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
		}

		$query_align = $trans_query_align;
		$match_align = $trans_match_align;

		# Check length threshold is met
		if (length($query_align) >= $match_length_threshold) {

#			# Check that the match hasn't already been extracted in reverse order
#			if (!exists($queries{"$match_name-$query_name"})) {
#
#				# Remove nucleotides to make length a multiple of 3
#				die "WHUT?\n" if (length($query_align) % 3 != 0);
#

				$count++;
	
				#my $name = "q_$query_name"."_t_$match_name";
			#	my $pair = {'QUERY_ALIGN' => $query_align,
			#				'MATCH_ALIGN' => $match_align,
			#				'LENGTH' => length($query_align)};
				my $pair = {'QUERY_NAME'  => $query_name,
				            'MATCH_NAME'  => $match_name,
					        'QUERY_ALIGN' => $query_align,
							'MATCH_ALIGN' => $match_align,
							'LENGTH' => length($query_align)};

				#$matches{$name} = $pair;
				#my %hit = ($name => $pair);

			#	my $should_store = int(rand(2));

			#	#undef(%hit);
				#%hit = ($$ => $count);
				until(Storable::is_storing() == 0){};
				store_fd($pair, $handle);

				#store_fd(\%hit, $handle) if $should_store == 1;
			#	print "$$ hit stored ($count)\n";
#
#				# Check if there is already a match between these two sequences
#				# if there is a match, the longer length one will be output
#				if (exists($matches{$name})) {
#					my $current_length = $matches{$name}->{'LENGTH'};
#					if ($current_length <= length($query_align)) {
#						$matches{$name} = $pair;
#					}
#				}
#				else {
#					$matches{$name} = $pair;
#				}
#				$queries{"$query_name-$match_name"}++;
#			}
		}
	}

	#return \%matches;
}

sub parse_self_blat_output {

	# Load transcriptome
	my %align = parse_fasta("$transcriptome.transdecoder.mRNA");

	# Check if we've already parsed blat output
	if (!-e $paralogs_seqs) {

		logger("Parsing self-blat output...");

		chomp(my $total_hits = `wc -l $blat_out`);
		$total_hits =~ s/(^\d+).*/$1/;

		my $lines_per_thread = ceil($total_hits / $max_forks);

		print "lines per file: $lines_per_thread\n";

		system("split $blat_out $split_blat_out_dir/$blat_out. -l $lines_per_thread");

		my @pids;
		my $select = new IO::Select;
		foreach my $file (glob("$split_blat_out_dir/*")) {

			my $pipe = new IO::Pipe;	
			my $pid = fork();

			# Child
			if ($pid == 0) {
				my $to_parent = $pipe->writer();
				#$to_parent->blocking(0);
				#$to_parent->blocking(1);
				#$to_parent->autoflush(1);

				#exit(0) if ($file !~ /\.aa$/ || $file !~ /\.ab/);
				#exit(0) if ($file !~ /\.aa$/);

				# Parse the results
				#my $hits = parse_self_blat_output_child($file, \%align);
				#my $hits = parse_self_blat_output_child($file, $to_parent, \%align);
				#parse_self_blat_output_child($file, $to_parent, \%align);
				parse_self_blat_output_child($file, $to_parent, \%align);

				#store_fd($hits, $to_parent);

			#	my %hash = (1 => $$);
			#	store_fd(\%hash, $to_parent);
				until(Storable::is_storing() == 0){};
				#my %hash = (DONE => $$);
				#store_fd(\%hash, $to_parent);
				store_fd({DONE => $$}, $to_parent);

				exit(0);
			}
			else {
				my $from_child = $pipe->reader();
				$select->add($from_child);
				push(@pids, $pid);
			}
		}
		#sleep(5);

		my $count = 0;

		my $id = 0;
		my %queries;
		my %matches;
		while (my @handles = $select->can_read) {

			my $handle = shift(@handles);
			my %hit = %{fd_retrieve($handle)};
			#use Data::Dumper;$Data::Dumper::Useqq=1;print Dumper \%hit;
			$count++;
			if (exists($hit{"DONE"})) {
				$select->remove($handle);
				$handle->close();
				next;
			}
			else {
				my $query_name = $hit{'QUERY_NAME'};
				my $match_name = $hit{'MATCH_NAME'};

				if (!exists($queries{"$match_name-$query_name"})) {

					my $name = "q_$query_name"."_t_$match_name";

					if (exists($matches{$name})) {
						my $current_length = $matches{$name}->{'LENGTH'};
						#if ($current_length <= length($query_align)) {
						if ($current_length <= $hit{'LENGTH'}) {
							#$matches{$name} = $pair;
							$matches{$name} = \%hit;
						}
					}
					else {
						#$matches{$name} = $pair;
						$matches{$name} = \%hit;
					}
					$queries{"$query_name-$match_name"}++;
				}
			}
		}
		print "$count hashes received\n";

		foreach my $pid (@pids) {
			waitpid($pid, 0);
		}

#		my $id = 0;
#		my %queries;
#		my %matches;
#		open(my $blat_out, "<", $blat_out);
#		while (my $line = <$blat_out>) {
#			chomp($line);
#
#			# Split the line on tabs and extract the information we want
#			my @line = split("\t", $line);
#			my ($query_name, $match_name, $query_align, $match_align) = 
#				($line[9], $line[13], $line[21], $line[22]);
#
#			# Don't want self matches
#			next if ($query_name eq $match_name);
#			
#			# Reverse translate amino acids to their corresponding nucleotides
#
#			my @query_align =  split(",", $query_align);
#			my @match_align =  split(",", $match_align);
#
#			my $query_nuc_align = $align{$query_name};
#			my $match_nuc_align = $align{$match_name};
#
#			my $trans_query_align;
#			foreach my $align (@query_align) {
#				$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
#			}
#
#			my $trans_match_align;
#			foreach my $align (@match_align) {
#				$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
#			}
#
#			$query_align = $trans_query_align;
#			$match_align = $trans_match_align;
#
#			# Check length threshold is met
#			if (length($query_align) >= $match_length_threshold) {
#
#				# Check that the match hasn't already been extracted
#				if (!exists($queries{"$match_name-$query_name"})) {
#
#					# Remove nucleotides to make length a multiple of 3
#					die "WHUT?\n" if (length($query_align) % 3 != 0);
#
#					my $name = "q_$query_name"."_t_$match_name";
#					my $pair = {'QUERY_ALIGN' => $query_align,
#								'MATCH_ALIGN' => $match_align,
#								'LENGTH' => length($query_align)};
#
#					# Check if there is already a match between these two sequences
#					# if there is a match, the longer length one will be output
#					if (exists($matches{$name})) {
#						my $current_length = $matches{$name}->{'LENGTH'};
#						if ($current_length <= length($query_align)) {
#							$matches{$name} = $pair;
#						}
#					}
#					else {
#						$matches{$name} = $pair;
#					}
#					$queries{"$query_name-$match_name"}++;
#				}
#			}
#		}

		# Prepare for output
		my @output;
		foreach my $key (sort { $a cmp $b} keys %matches) {
			my $pair = $matches{$key};
			push(@output, ">$key\n");
			push(@output, "$pair->{QUERY_ALIGN}\n");
			push(@output, "$pair->{MATCH_ALIGN}\n\n");
			$id++;
		}

		if ($id == 0) {
			logger("No blat hits met the requirements.\n");
			exit(0);
		}
		logger("Completed parsing blat output, $id blat hit(s) met the requirements.\n");

		open(my $output_file, ">", $paralogs_seqs);
		print {$output_file} @output;
		close($output_file);
		die;
	}
	else {
		logger("Using previously parsed blat output in '$paralogs_seqs'.");
	}
}

sub run_kaks_calc {

	# Check if we've run KaKs_Calculator before
	if (!-e $paralogs_ka_ks_out) {
		# Run KaKs_Calculator to get Ks values
		logger("Running KaKs_Calculator...");
		system("$kaks_calculator -i $paralogs_seqs -o $paralogs_ka_ks_out -m $model >/dev/null") && die;
		logger("Completed KaKs_Calculator.\n");
	}
	else {
		logger("Using KaKs_Calculator output from previous run in '$paralogs_ka_ks_out'.");
	}
}

sub parse_ks_values {

	# Check if we've parsed Ks values before
	if (!-e $final_ks_values) {

		# Open KaKs_Calculator output, parse Ks values and convert to .csv for R
		open(my $ka_ks_calculator_output, "<", $paralogs_ka_ks_out);
		open(my $ks_csv, ">", $final_ks_values);
		while (my $line = <$ka_ks_calculator_output>) {
			chomp($line);
			my @line = split(/\s+/, $line);

			# Skip first line containing headers
			if ($. == 1) {
				print {$ks_csv} "ks\n";
				next;
			}

			my $ks = $line[3];
			$ks = 0 if ($ks eq "NA");
			#next if ($ks == 0 && $exclude_zero);

			print {$ks_csv} "$ks\n"
		}
		close($ks_csv);
		close($ka_ks_calculator_output);
	}
}

sub create_ks_plot {

	# Create PDF plot of output
	logger("Creating Ks plot in R...");
	if ($exclude_zero) {
		if ($bin_size == 0) {
			system("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values'); 
				data <- data\\\$ks[data\\\$ks <= $ks_max & data\\\$ks > 0]; 
				plot(density(data), main=expression(paste('K'[s], ' Density Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
		else {
			system("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values'); 
				data <- data\\\$ks[data\\\$ks <= $ks_max & data\\\$ks > 0]; 
				hist(data, breaks=seq($ks_min,$ks_max,by=$bin_size), 
					main=expression(paste('K'[s], ' Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
	}
	else {
		if ($bin_size == 0) {
			system("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values'); 
				data <- data\\\$ks[data\\\$ks <= $ks_max]; 
				plot(density(data), main=expression(paste('K'[s], ' Density Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
		else {
			system("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values'); 
				data <- data\\\$ks[data\\\$ks <= $ks_max]; 
				hist(data, breaks=seq($ks_min,$ks_max,by=$bin_size), 
					main=expression(paste('K'[s], ' Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
	}
	logger("\nCompleted Ks plot.\n");
}

sub reverse_translate {
	my $settings = shift;

	my $dna = $settings->{'DNA'};
	my $prot = $settings->{'PROT'};

	my $regex;
	foreach my $index (0 .. length($prot) - 1) {
		my $char = substr($prot, $index, 1);
		$regex .= $rev_codon_table{$char};
	}

	my $translation;
	if ($dna =~ /($regex)/) {
		$translation = $1;
	}
	elsif (reverse($dna) =~ /($regex)/) {
		$translation = $1;
	}
	else {
		die "Protein sequence could not be reverse translated.\n";
	}

	return $translation;
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec);
	}

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

sub parse_fasta {
	my $filename = shift;

	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";
	chomp(my @data = <$alignment_file>);
	close($alignment_file);
	
	my $taxon;
	foreach my $line (@data) {
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			$taxon =~ s/-/_/g;
			$align{$taxon} .= $line;
		}
	}
	return %align;
}

sub logger {
	my $msg = shift;

	my $time = "[".localtime(time())."]";

	# Allow for new lines and carriage returns before message
	if ($msg =~ s/^\n//) {
		$time = "\n$time";
	}
	elsif ($msg =~ s/^\r//) {
		$time = "\r$time";
	}

	print "$time $msg\n"; 
}

sub INT_handler {
	kill(15, $pid);
	exit(0);
}

sub TERM_handler {

	# We only want to clean up if the initial input was a sequence file
	if (!$input_is_dir) {
		logger("\rKeyboard interrupt detected, stopping analyses and cleaning up.");

		# Move back into the directory script was called in
		chdir("..");

		# Try to delete directory five times, if it can't be deleted print an error message
		my $count = 0;
		until (!-e $project_name || $count == 5) {
			$count++;

			remove_tree($project_name, {error => \my $err});
			sleep(1);
		}
		logger("Could not clean all files in './$project_name/'.") if ($count == 5);
	}
	else {
		logger("\rKeyboard interrupt detected, stopping analyses.");
	}

	exit(0);
}

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}

sub usage {
	return "Usage: perl plot-ks.pl [TRANSCRIPTOME] [OPTIONS]...\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Generate a Ks plot for a given transcriptome in fasta format

  -m, --model                       model used by KaKs_Calculator to determine Ks (default: YN)
  -t, --match_length_threshold      the minimum number of basepairs the matching sequences must be (default: 300 bp)
  -x, --exclude_zero                used to exclude Ks = 0 from plot, useful for Trinity transcriptomes
  -b, --bin_size                    size of bins used in histogram of Ks plot, set to 0 for a density plot (default: 0.05)
  --ks_min                          lower boundary for x-axis of Ks plot (default: Ks = 0)
  --ks_max                          upper boundary for x-axis of Ks plot (default: Ks = 3)
  -h, --help                        display this help and exit

Examples:
  perl plot-ks.pl assembly.fa -b 0.01               generates a Ks (YN model) plot from [0, 3] using a bin size of 0.01, and contigs 
                                                    with at least 300 bp of homologous sequence

  perl plot-ks.pl assembly.fa -x -m NG              generates a Ks (NG model) plot from (0, 3] using a bin size of 0.05, and contigs 
                                                    with at least 300 bp of homologous sequence

  perl plot-ks.pl assembly.fa --ks_max 5 -t 500     generates a Ks (YN model) plot from [0, 5] using a bin size of 0.01, and contigs 
                                                    with at least 500 bp of homologous sequence

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
