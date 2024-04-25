use anyhow::{anyhow, Context, Result};
use camino::Utf8PathBuf as PathBuf;
use clap::{App, Arg};
use csv::Reader;
use csv::Writer;
use futures::future::join_all;
use regex::Regex;
use reqwest::Client;
use std::fs::{create_dir_all, File};
use std::{collections::HashSet, fs};

async fn fetch_full_filename(client: &Client, accession: &str) -> Result<(String, String)> {
    let (db, acc) = accession
        .trim()
        .split_once('_')
        .ok_or_else(|| anyhow!("Invalid accession format"))?;
    let (number, _) = acc.split_once('.').unwrap_or((acc, "1"));
    let number_path = number
        .chars()
        .collect::<Vec<_>>()
        .chunks(3)
        .map(|chunk| chunk.iter().collect::<String>())
        .collect::<Vec<_>>()
        .join("/");

    let base_url = format!(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}",
        db, number_path
    );
    let directory_response = client.get(&base_url).send().await?;
    if !directory_response.status().is_success() {
        return Err(anyhow!(
            "Failed to open genome directory: HTTP {}",
            directory_response.status()
        ));
    }
    let text = directory_response.text().await?;
    let link_regex = Regex::new(r#"<a href="([^"]*)""#)?;

    for cap in link_regex.captures_iter(&text) {
        let name = &cap[1];
        // Check if name ends with '/', and remove it if so
        let clean_name = if name.ends_with('/') {
            name.strip_suffix('/').unwrap()
        } else {
            name
        };
        if clean_name.starts_with(db)
            && clean_name
                .split('_')
                .nth(1)
                .map_or(false, |x| x.starts_with(number))
        {
            return Ok((format!("{}/{}", base_url, clean_name), clean_name.into()));
        }
    }

    Err(anyhow!(
        "No matching genome found for accession {}",
        accession
    ))
}

async fn download_with_retry(
    client: &Client,
    url: &str,
    file_name: PathBuf,
    retry_count: u32,
) -> Result<()> {
    let mut attempts = retry_count;
    while attempts > 0 {
        let response = client.get(url).send().await;
        match response {
            Ok(resp) if resp.status().is_success() => {
                let data = resp
                    .bytes()
                    .await
                    .context("Failed to read bytes from response")?;
                fs::write(file_name, &data).context("Failed to write data to file")?;
                return Ok(());
            }
            _ => {
                eprintln!("Failed to download file: {}. Retrying...", url);
                attempts -= 1;
            }
        }
    }

    Err(anyhow!(
        "Failed to download file after {} retries: {}",
        retry_count,
        url
    ))
}

async fn process_accession(
    client: &Client,
    accession: &str,
    location: &PathBuf,
    retry: Option<u32>,
) -> Result<()> {
    let retry_count = retry.unwrap_or(3); // Default retry count

    let (base_url, full_name) = fetch_full_filename(client, accession).await?;

    let suffixes = vec!["_genomic.fna.gz", "_protein.faa.gz"]; //, "_assembly_report.txt"];
    let standalone = vec!["md5checksums.txt"];

    for suffix in suffixes.iter() {
        let url = format!("{}/{}{}", base_url, full_name, suffix); // Correctly format the URL for each file type
        let file_name = format!("{}{}", accession, suffix); // Generate file name using the directory name and suffix
        let path = location.join(&file_name); // Create the full path for the file
        download_with_retry(client, &url, path, retry_count).await?;
    }

    // download standalone files (mostly md5checksums.txt)
    for filename in standalone {
        let url = format!("{}/{}", base_url, filename);
        let file_name = format!("{}_{}", accession, filename); // Generate file name using the directory name and suffix
        let path = location.join(&file_name); // Create the full path for the file
        download_with_retry(client, &url, path, retry_count).await?;
    }

    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    let matches = App::new("assembly-downloader")
        .version("1.0")
        .author("Tessa Pierce-Ward <ntpierce@gmail.com>")
        .about("Download NCBI Assembly Datasets from Accession List")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .takes_value(true)
                .help("Input CSV file containing accession numbers in the first column")
                .required(true),
        )
        .arg(
            Arg::new("failed")
                .short('f')
                .long("failed")
                .takes_value(true)
                .help("Output file to write failed accession numbers")
                .required(true),
        )
        .arg(
            Arg::new("retry-times")
                .short('r')
                .long("retry-times")
                .takes_value(true)
                .default_value("3")
                .help("Number of times to retry a failed download"),
        )
        .arg(
            Arg::new("location")
                .short('l')
                .long("location")
                .takes_value(true)
                .default_value(".")
                .help("Directory location where files will be downloaded"),
        )
        .get_matches();

    let input_csv = matches.value_of("input").unwrap();
    let failed_csv = matches.value_of("failed").unwrap();
    let retry_times: u32 = matches.value_of_t("retry-times")?; // Safely parsing the retry times as u32
    let location = matches.value_of("location").unwrap();

    let download_path = PathBuf::from(location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }

    // Open the file containing the accessions synchronously
    let file = File::open(input_csv)?;
    let mut rdr = Reader::from_reader(file);

    // Initialize a HashSet to store unique accessions
    let mut accessions = HashSet::new();

    // Read all accessions into the HashSet to remove duplicates
    for result in rdr.deserialize::<(String,)>() {
        let record = result?;
        let accession = record.0.trim();
        if !accession.is_empty() {
            accessions.insert(accession.to_string());
        }
    }

    // Initialize the HTTP client
    let client = Client::new();
    let mut failed_writer = Writer::from_path(failed_csv)?;
    failed_writer.write_record(&["accession", "url"])?;

    // Collect accessions into a vector for easier chunking
    let accessions_vec: Vec<_> = accessions.iter().collect();

    // NCBI rate-limits to 3 requests/second, so process 3 at a time:
    for chunk in accessions_vec.chunks(3) {
        let futures = chunk.iter().map(|accession| {
            let client_ref = &client;
            let download_path_ref = &download_path;
            let accession_clone = accession.to_owned();
            async move {
                match process_accession(
                    client_ref,
                    &accession_clone,
                    download_path_ref,
                    Some(retry_times),
                )
                .await
                {
                    Ok(_) => Ok(accession_clone),
                    Err(e) => Err((accession_clone, e)),
                }
            }
        });

        // Wait for all accessions in the current chunk to be processed
        let results = join_all(futures).await;

        for result in results {
            match result {
                Ok(accession) => println!("Successfully processed accession: {}", accession),
                Err((accession, e)) => {
                    let err_message = e.to_string();
                    let parts: Vec<&str> = err_message.split("retries: ").collect();
                    let failed_url = parts.get(1).unwrap_or(&"Unknown URL").trim();

                    failed_writer.write_record(&[accession, failed_url])?;
                    eprintln!(
                        "Failed to process accession: {}. Error: {}",
                        accession, err_message
                    );
                }
            }
        }

    }

    Ok(())
}
