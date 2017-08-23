#use "topfind";;
#thread;;
#require "ketrew";;
#require "coclobas.ketrew_backend";;
open Printf;;

let () = Unix.putenv "KETREW_CONFIGURATION" (Sys.getenv "KETREW_ROOT" ^ "/configuration.ml");;

let pythonpath = String.concat ":" [
  "/tmp/cartpipeline/src/pyhpa";
  "/tmp/cartpipeline/src/pycgds";
  "/tmp/cartpipeline/src/pyagg";
];;
let projectpath = "/Users/eczech/projects/hammer/";;

let tcga_studies =
  let parseargs argv = match argv with
    | [] | [_] -> raise (Invalid_argument "No arguments given (must pass TCGA study id names)")
    | hd::tl -> tl
  in parseargs (Array.to_list Sys.argv);;

print_string ("Running pipeline with TCGA study ids: " ^ (String.concat ", " tcga_studies));;

let run_pipeline cmd1 =
  let open Ketrew.EDSL in
  let prjpath path =  projectpath ^ path in
  let volume_mounts = [
    `Local (prjpath "repos/cartpipeline/python", "/tmp/cartpipeline/src");
    `Local (prjpath "cache/pipeline", "/tmp/cartpipeline/data");
    `Local (prjpath "cache/tcga", "/tmp/tcgacache")
  ] in
  let dockerize pgm =
    Coclobas_ketrew_backend.Plugin.local_docker_program
      ~tmp_dir: "/tmp/secotrec-local-shared-temp"
      ~base_url: "http://coclo:8082"
      ~image: "py3.5-v1"
      ~volume_mounts: volume_mounts
      pgm in
  let pypgm script args =
    Program.(
      shf "export PYTHONPATH=\"%s:$PYTHONPATH\"" pythonpath &&
      shf "python /tmp/cartpipeline/src/%s %s" script args
    ) in
  let collect_meta =
    workflow_node without_product
      ~name:"Gene Meta Generation"
      ~make: (dockerize (pypgm
        "pyhpa/script/gene_selector.py"
        "--output /tmp/cartpipeline/data/gene_meta.csv"
      )) in
  let collect_expression = List.map (fun study_id ->
    workflow_node without_product
      ~name:(sprintf "Expression Data Collection (%s)" study_id)
      ~edges:[depends_on collect_meta]
      ~make: (dockerize (pypgm
        "pycgds/script/tcga_expression.py" (sprintf
          "--study-id %s \
          --output /tmp/cartpipeline/data/expression_data_%s.csv \
          --gene-meta-path /tmp/cartpipeline/data/gene_meta.csv \
          --cache-dir /tmp/tcgacache" study_id study_id
        )
      ))
    ) tcga_studies in
  let expression_paths = String.concat " " (
    List.map (fun study_id ->
      sprintf "/tmp/cartpipeline/data/expression_data_%s.csv" study_id
    ) tcga_studies
  ) in
  let aggregate =
    workflow_node without_product
      ~name:"Result Aggregation"
      ~edges:(List.map depends_on collect_expression)
      ~make: (dockerize (pypgm
        "pyagg/script/aggregation.py" (sprintf
          "--gene-meta-path /tmp/cartpipeline/data/gene_meta.csv \
          --gene-exp-paths %s \
          --output /tmp/cartpipeline/data/pipeline_result.csv"
          expression_paths
        )
      ))
  in Ketrew.Client.submit_workflow aggregate;;

let () = run_pipeline "ls";;

(* utop cart_pipeline_v3.ml brca_tcga prad_tcga luad_tcga skcm_tcga laml_tcga gbm_tcga lgg_tcga coadread_tcga paad_tcga ov_tcga kirc_tcga meso_tcga *)
