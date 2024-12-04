# Changelog

## [1.4.0](https://github.com/3d-omics/mt_quant/compare/v1.3.0...v1.4.0) (2024-12-04)


### Features

* add slurm, cache, test scripts ([4ab571f](https://github.com/3d-omics/mt_quant/commit/4ab571feea25d49cdc85dcc1be2f0b98b722e30f))
* bracken ([8183326](https://github.com/3d-omics/mt_quant/commit/8183326760c5a1cc6e605d6482a7d10883f18763))
* cache rule ([e7ac7bd](https://github.com/3d-omics/mt_quant/commit/e7ac7bd211bf1806191ca14282f02c5934db050e))
* configfile ([ebac6c0](https://github.com/3d-omics/mt_quant/commit/ebac6c0c678307df2f506af7bc9be0e550ee14dd))
* do bracken on all taxonomic levels, only show species in multiqc ([d8efee3](https://github.com/3d-omics/mt_quant/commit/d8efee335a5ba8a613eef7b741dcf4f4ac80b009))
* get the aggregated counts directly from subread ([7f8900e](https://github.com/3d-omics/mt_quant/commit/7f8900e281dc46356f46270f71b2476be7b89a47))
* use col1 and col2 to see if it reduces ram usage ([2f527c2](https://github.com/3d-omics/mt_quant/commit/2f527c23406c0facc472e6b0ca2b9c4776f29fe0))
* use subread from wrappers ([e9971d8](https://github.com/3d-omics/mt_quant/commit/e9971d830237348498a0746536168dc81c90f8e5))


### Bug Fixes

* add all the bracken files to the copy in /dev/shm ([c4df691](https://github.com/3d-omics/mt_quant/commit/c4df6917a798b304b4d7fd8c5a916b08fcd896a8))
* dont use slurm within slurm , dumbass ([acf4c57](https://github.com/3d-omics/mt_quant/commit/acf4c579fb3638f49c82be70866d2ecf54d10440))
* get the correct featurecounts file ([a2beada](https://github.com/3d-omics/mt_quant/commit/a2beada0c3785b66b720cce0c89456ecc2c0df55))
* handle manually subread sample renaming; corner case in csvtk ([0a15ae8](https://github.com/3d-omics/mt_quant/commit/0a15ae848f50e1bf3fb198bd34830d23cb497c59))
* put clean in the preprocess group ([24be3b8](https://github.com/3d-omics/mt_quant/commit/24be3b86f13b71422aad07a8f411e243bfff212d))
* put fastq in the correct group ([e6a3462](https://github.com/3d-omics/mt_quant/commit/e6a34626d6ad9cab0794af9ab4960fcd29e4d728))
* put slurm executor ([f4b65f6](https://github.com/3d-omics/mt_quant/commit/f4b65f66fec74287ed559ea5997f1b2088d81da6))
* put tmp folders for index and map ([a3c11ab](https://github.com/3d-omics/mt_quant/commit/a3c11abbb3efb85bd7d73211a6b8275a10703d8c))
* remove more mentions to coverm ([900e667](https://github.com/3d-omics/mt_quant/commit/900e6678e3cee404195ff3c2d2c27f1ef849d2b8))
* remove tmp dirs if they exist ([ccf4e09](https://github.com/3d-omics/mt_quant/commit/ccf4e094bf3c8f2fc8b7821edd21700b41249904))
* rerun incomplete ([99941d1](https://github.com/3d-omics/mt_quant/commit/99941d10d326b67bff1e58ab7ff01251ba70066f))
* undo bracken copy ([47f72eb](https://github.com/3d-omics/mt_quant/commit/47f72ebc463093037ff202db18801c375fb68d9a))
* use proper featurecounts file ([d76c680](https://github.com/3d-omics/mt_quant/commit/d76c68058a50814728863d866cc7678cd4b4d75e))


### Performance Improvements

* cache but omit software versions ([09eb311](https://github.com/3d-omics/mt_quant/commit/09eb3112d874f7a09a7101a62d5eefed3e8c8174))
* disable groups ([6d17ce8](https://github.com/3d-omics/mt_quant/commit/6d17ce88c3a7c277b016a51b17084aead5daf505))
* double ram on subread join ([f3c356f](https://github.com/3d-omics/mt_quant/commit/f3c356f315284489e64020ac5d4bf7af6d47b116))
* put reads in a group and disable localrule ([e7c2ad7](https://github.com/3d-omics/mt_quant/commit/e7c2ad7d3d6fa1195cbd4e455f66a91c0b37bdf6))
* raise ram to tsv join ([a06dcad](https://github.com/3d-omics/mt_quant/commit/a06dcad6128df5bd0ca7fb8dd55d3023de731fae))
* **subread:** raise ram to 4GB ([b0ea99a](https://github.com/3d-omics/mt_quant/commit/b0ea99a72b2e2264e4fdd3fe63f5ffe29abea477))

## [1.3.0](https://github.com/3d-omics/mt_quant/compare/v1.2.0...v1.3.0) (2024-11-29)


### Features

* **coverm:** mark individual files as temporary ([25e4c86](https://github.com/3d-omics/mt_quant/commit/25e4c8625e0d0eece98c16fb382356324edf01ac))
* **coverm:** rearrange files again ([8058735](https://github.com/3d-omics/mt_quant/commit/8058735f11e9eac0cf20cc6eb66c30aff5c011d2))
* **coverm:** remove functions ([015c4f8](https://github.com/3d-omics/mt_quant/commit/015c4f8a76254ab25219269ea50ba7e64b5aaff4))
* **htseq:** do not compress output and make it temp ([6d9159c](https://github.com/3d-omics/mt_quant/commit/6d9159c12176d656edc03743ff4f66bbdddcc956))
* remove r packages, update and pin environments ([5785f85](https://github.com/3d-omics/mt_quant/commit/5785f859524628689b1ae974c03f849e9743d738))
* restructure bowtie and use csvtk to join files ([acc9ddf](https://github.com/3d-omics/mt_quant/commit/acc9ddfbf35704d709d67ec2cb431f026215f5b6))
* **subread:** introduce middle step ([365d41e](https://github.com/3d-omics/mt_quant/commit/365d41e479b2ab552755e7dd12dcff410b1da50e))
* **subread:** quantify trnas and rrnas too ([d459a74](https://github.com/3d-omics/mt_quant/commit/d459a7417ca304181383b4606909411a50349fb5))
* **subread:** solve dependency; add reports to multiqc ([2b7b9a3](https://github.com/3d-omics/mt_quant/commit/2b7b9a385e4a34ff09a74b77cbf17085921f7374))
* substitute custom r script for csvtk ([0cacb96](https://github.com/3d-omics/mt_quant/commit/0cacb9600d8ad4fcf4fba36b7461a9ed325303fa))


### Bug Fixes

* **bowtie2:** small index ([7f17e44](https://github.com/3d-omics/mt_quant/commit/7f17e447dd149e470d6659cd44e6bcb9c1ee1e9e))
* give multiqc the correct files ([1f35182](https://github.com/3d-omics/mt_quant/commit/1f35182fdb440fe8b831ef3a2722c56deb276cdc))
* **htseq:** remove mention to file ([53465d9](https://github.com/3d-omics/mt_quant/commit/53465d911bddcb66ccdde781e0246d5cbbf0f395))
* remodel to take the gff from dram ([7a5d362](https://github.com/3d-omics/mt_quant/commit/7a5d362870f9f713473f6408789ab855e861e0c3))

## [1.2.0](https://github.com/3d-omics/mt_quant/compare/v1.1.2...v1.2.0) (2024-11-12)


### Features

* coverm++ ([4b170ac](https://github.com/3d-omics/mt_quant/commit/4b170ac9d0053b08a2effae8666fde34cd78f7c1))
* demand snakemake 8 ([5cec580](https://github.com/3d-omics/mt_quant/commit/5cec58083427c5154cd995b8f0a0254d6d885238))
* disable crams in prerpocessing ([f329758](https://github.com/3d-omics/mt_quant/commit/f329758bd7520a1d1be4e2d9fbd6d1a168684256))
* disable crams in quantify ([c618038](https://github.com/3d-omics/mt_quant/commit/c61803828b8730bbdc8c4037293cfe104491f1bb))
* group jobs by sample, update profile ([658ae38](https://github.com/3d-omics/mt_quant/commit/658ae386e3a9859f650d8cc590419018e8eeea75))
* import helpers. insert reads into prerpocess. remove per library report ([af3f230](https://github.com/3d-omics/mt_quant/commit/af3f230292480d39c56b9250e60722686b6a9a9a))
* insert hosts into preprocess ([4c285c0](https://github.com/3d-omics/mt_quant/commit/4c285c0d51aa621273adb67d9c5a790511a9c6e2))
* multqc handles more files ([4d44cad](https://github.com/3d-omics/mt_quant/commit/4d44cad32c5260c2b055a10f575bbe99e8c45216))
* point to base env and clean other ([b1b34e0](https://github.com/3d-omics/mt_quant/commit/b1b34e07b33f9d414687428f94d17e6c2b9ccb50))
* refactor fastp. use wrappers ([da1a88d](https://github.com/3d-omics/mt_quant/commit/da1a88d286fb01ad024e802fbc8468eb367d442f))
* refactor quantify ([e63eca6](https://github.com/3d-omics/mt_quant/commit/e63eca681e4017f41b729bd615371e350479b1e8))
* refactor quantify ([d59ce09](https://github.com/3d-omics/mt_quant/commit/d59ce099c04918e617ab238f78d64a9730d54b3d))
* refactor ribodetector ([cee57ef](https://github.com/3d-omics/mt_quant/commit/cee57ef4f58868b568384cd736166091addb4605))
* refactor ribodetector ([bc4f921](https://github.com/3d-omics/mt_quant/commit/bc4f921f6f9333e3fa742a00136ad50432615dd8))
* remove fastp fastqcs from report ([9323aec](https://github.com/3d-omics/mt_quant/commit/9323aeca487543bcc2cd30f5700533a4fb3238a5))
* remove fastp functions ([5eda525](https://github.com/3d-omics/mt_quant/commit/5eda52533d23fe86beee81f0c58ea25aa70e7f42))
* separate bowtie2 env from quant ([71a54a6](https://github.com/3d-omics/mt_quant/commit/71a54a634636a2e72fea9759b5914d3b5c07e2e4))
* update quantify rule names, fix multiqc output ([96cb3e8](https://github.com/3d-omics/mt_quant/commit/96cb3e8436d9bacf88751e2cd5839787c7a7c6fc))


### Bug Fixes

* **fastp:** correct group name ([2751414](https://github.com/3d-omics/mt_quant/commit/2751414963620108d653dfddde90bd8461192346))
* remove mention to ram ([e009781](https://github.com/3d-omics/mt_quant/commit/e009781995dc198d4e910f86fcf13035cf406f0d))
* **star:** point to u1, u2 ([f6344f1](https://github.com/3d-omics/mt_quant/commit/f6344f13b8bd4e807fa340d72d50ae2ec2f9343b))
* **star:** point to u1, u2 .fq ([fb68d67](https://github.com/3d-omics/mt_quant/commit/fb68d676317fcd2056f95958cd5ef1fc1006d4c5))


### Performance Improvements

* add fastq to group ([ee27676](https://github.com/3d-omics/mt_quant/commit/ee27676af4c3d44f385eb37e1064f9aa12e684e6))
* disable group in bowtie2 ([e3a26e1](https://github.com/3d-omics/mt_quant/commit/e3a26e1cbb911b0921f015cd775284399e2d6e3e))
* **fastp:** default extra string ([1c9dfc4](https://github.com/3d-omics/mt_quant/commit/1c9dfc417a04b806c023c8f33664da457840be99))
* **profile:** increase ram for star index ([06bc617](https://github.com/3d-omics/mt_quant/commit/06bc617046a5c868ac8127649c7f0c8eff9cd176))

## [1.1.2](https://github.com/3d-omics/mt_quant/compare/v1.1.1...v1.1.2) (2024-10-04)


### Bug Fixes

* decompress mag gtf ([db32655](https://github.com/3d-omics/mt_quant/commit/db326551766029df7006b34f6637a85f4e3e868e))

## [1.1.1](https://github.com/3d-omics/mt_quant/compare/v1.1.0...v1.1.1) (2024-10-01)


### Bug Fixes

* **preprocess/index:** make params.prefix str so caching works ([2a8a205](https://github.com/3d-omics/mt_quant/commit/2a8a205bb9f506c0716ab954bee1e90db7e3056d))

## [1.1.0](https://github.com/3d-omics/mt_quant/compare/v1.0.1...v1.1.0) (2024-08-06)


### Features

* pin dependencies ([1d39136](https://github.com/3d-omics/mt_quant/commit/1d3913637838651520d2a5f7c69820ea1e017b10))

## [1.0.1](https://github.com/3d-omics/mt_quant/compare/v1.0.0...v1.0.1) (2024-06-17)


### Performance Improvements

* raise ram for aggregators ([258ec0f](https://github.com/3d-omics/mt_quant/commit/258ec0f5e3066ffa680247cf6766239f5c3e4116))

## [1.0.0](https://github.com/3d-omics/mt_quant/compare/0.0.1...v1.0.0) (2024-06-14)


### ⚠ BREAKING CHANGES

* release 1.0.0

### Features

* add gha updater ([6dd2b43](https://github.com/3d-omics/mt_quant/commit/6dd2b43bb7d3ed4bd4f3c8942bf8a1f1abc12c90))
* add per gene counts with featureCount ([3f8fd05](https://github.com/3d-omics/mt_quant/commit/3f8fd05688d09444ba8ec70d9a28739df4b7c1e2))
* **bowtie2:** write index within the rule ([f171f0c](https://github.com/3d-omics/mt_quant/commit/f171f0cfca075181830cc1d87de92e966ac85dd3))
* **bowtie2:** write index within the rule ([5f17b3f](https://github.com/3d-omics/mt_quant/commit/5f17b3f3ee87f20b56b7c88fb07e76fea4357f9b))
* compress all count data ([72e9a9e](https://github.com/3d-omics/mt_quant/commit/72e9a9e9c4579c9739fcad053962e3587012343d))
* compute counts over genes ([a222296](https://github.com/3d-omics/mt_quant/commit/a2222960e94e70a3e699e3e3f9de6a24163f1da2))
* **coverm:** coverm + cram ([dbe1b09](https://github.com/3d-omics/mt_quant/commit/dbe1b09b9f013bcc915bf41799b5e695faee5535))
* **coverm:** coverm + cram ([869df57](https://github.com/3d-omics/mt_quant/commit/869df5702d65fc699c8a4f7b5fdd9cab3f28907d))
* **coverm:** delete dead code. -F 4 ([0661d62](https://github.com/3d-omics/mt_quant/commit/0661d62ca53ea777f63668b00061f5215214f776))
* **coverm:** delete dead code. -F 4 ([dea6a0a](https://github.com/3d-omics/mt_quant/commit/dea6a0a98869d59f7135c5061b946ef292b755be))
* **fastp:** remove process substitution and use level 9 compression ([7a685dd](https://github.com/3d-omics/mt_quant/commit/7a685ddb0e79a52f2f852c109de832e4ca3b007a))
* **fastp:** remove process substitution and use level 9 compression ([a45a22f](https://github.com/3d-omics/mt_quant/commit/a45a22f5ca29fac1fb19203f097a1718f86bc0ab))
* **fastp:** remove process substitution and use level 9 compression ([6834cc5](https://github.com/3d-omics/mt_quant/commit/6834cc5c2797e27af0aa62244933f6b5cd564da7))
* **fastp:** use process substituttion and pigz -11 ([9b956cb](https://github.com/3d-omics/mt_quant/commit/9b956cbbc0fe7964c24d2c068513eba1ed48dccc))
* **fastp:** use process substituttion and pigz -11 ([4f3b0b9](https://github.com/3d-omics/mt_quant/commit/4f3b0b98f43eb8f81140ad9dad73eb127674531f))
* **fastp:** use process substituttion and pigz -11 ([c49b538](https://github.com/3d-omics/mt_quant/commit/c49b5384099557ee6d989f20ee3422cd66b3c34b))
* **features:** add the usual suspects ([a2ec8c4](https://github.com/3d-omics/mt_quant/commit/a2ec8c4b6b2a21c1bee50006423c4f38da300ec9))
* **features:** move databases and mag catalogue. @ in mag mock ([f9c5f54](https://github.com/3d-omics/mt_quant/commit/f9c5f549a8c34ed35de28db30de68ac8fd45eea1))
* get gene quantifications with htseq ([1a461e0](https://github.com/3d-omics/mt_quant/commit/1a461e0879f42ee95bb086f08c4f33554a80e29d))
* ignore rstudio stuff ([651d20e](https://github.com/3d-omics/mt_quant/commit/651d20edf85eb4ddfa3f83a84f90989558aa1a52))
* ignore rstudio stuff ([39767ac](https://github.com/3d-omics/mt_quant/commit/39767ac343d11158083f946f0353618ea14349d6))
* **params:** add @ as the default mag-contig separator ([e33b82d](https://github.com/3d-omics/mt_quant/commit/e33b82dd02fb62170fc202ec1721e9762209cb35))
* **params:** add common coverm params ([47f0caf](https://github.com/3d-omics/mt_quant/commit/47f0cafd143982239d511d35f9bbcd4f9e41e316))
* **params:** add common coverm params ([7a00763](https://github.com/3d-omics/mt_quant/commit/7a007635e6897a4f15d88215fa48426cc8bde8a0))
* **pre-commit:** add conventional-pre-commit ([b6df867](https://github.com/3d-omics/mt_quant/commit/b6df8674d963a2ba7d1fe0583828194df1624ecb))
* **pre-commit:** add R linters ([b31ddec](https://github.com/3d-omics/mt_quant/commit/b31ddecda5e91fc55900cbcbe8b41ef2b98a9e3b))
* **pre-commit:** autoupgrade ([99ec679](https://github.com/3d-omics/mt_quant/commit/99ec679e711dc425bf487f14e4c79868ce31d72a))
* **pre-commit:** autoupgrade ([e156263](https://github.com/3d-omics/mt_quant/commit/e156263a69667677f4760bf018590d05728762ce))
* **pre-commit:** update all hooks ([e7cf184](https://github.com/3d-omics/mt_quant/commit/e7cf18403a6306b994d32b35ddf61344ad3583e1))
* **preprocess:** set adaptors in samples. refactor. Update readme ([649d019](https://github.com/3d-omics/mt_quant/commit/649d019b6e6f6ce8e7bfb789767a43c213876969))
* **quantify:** refactor ([be70cab](https://github.com/3d-omics/mt_quant/commit/be70cabfb377b6b8237ead3d4726ae74d1ad3172))
* **quantify:** refactor coverm functions ([6a33bc0](https://github.com/3d-omics/mt_quant/commit/6a33bc069070d80a994bda0076ec107bbe41c7c7))
* **README:** update repo name ([d516776](https://github.com/3d-omics/mt_quant/commit/d516776f75ac68061385800ffbc22120942eab08))
* **README:** update repo name ([c62c2be](https://github.com/3d-omics/mt_quant/commit/c62c2bead2af9c64f7388df5add11fc581db2ab6))
* **reads:** refactor file and adaptor functions ([ba42fec](https://github.com/3d-omics/mt_quant/commit/ba42fec6e1849b4317129ab10ecfecde04cdbea2))
* **reads:** split into subworkflows ([d6bc3ba](https://github.com/3d-omics/mt_quant/commit/d6bc3ba79d5c5f68ea988b2eeed38041f4056ed8))
* **reads:** split into subworkflows ([1616a6f](https://github.com/3d-omics/mt_quant/commit/1616a6fee0eac01e1cdd3327aea7eefd5f9f738c))
* **report:** add the definitive report rule ([a708af3](https://github.com/3d-omics/mt_quant/commit/a708af34180038491648991011c8181eae133ea1))
* **report:** raise ram ([f1aaf9d](https://github.com/3d-omics/mt_quant/commit/f1aaf9dc882f918631335b0f6e08d21c24f83836))
* **report:** raise ram ([ae130be](https://github.com/3d-omics/mt_quant/commit/ae130bef10146a505e157829a44d52840fa09fae))
* **ribodetector:** remove temp ([86c5ac2](https://github.com/3d-omics/mt_quant/commit/86c5ac2a432bd31fdc04f44af833b6de9ff762b9))
* **ribodetector:** remove temp ([0617fd3](https://github.com/3d-omics/mt_quant/commit/0617fd3b0d217a7a4a4937abf3929ee2482f8b41))
* **star:** compress within the rule, raise compression level ([eca2288](https://github.com/3d-omics/mt_quant/commit/eca2288bf7c5d7b85b2a7c3ac4fed052d3bb661a))
* **star:** compress within the rule, raise compression level ([d93087b](https://github.com/3d-omics/mt_quant/commit/d93087b4743aee0075380f58f3fe77e870ad33b8))
* **star:** refactor functions ([6258c9b](https://github.com/3d-omics/mt_quant/commit/6258c9bc70b8f69db39e539c4951d5585d656414))
* **star:** refactor functions ([85b0467](https://github.com/3d-omics/mt_quant/commit/85b0467e5560ed1c72e46b149866c6285d3475a2))
* the big rename ([0d9833f](https://github.com/3d-omics/mt_quant/commit/0d9833fe902111347c15a5afe802426511a09fb0))
* update github actions ([2e9d284](https://github.com/3d-omics/mt_quant/commit/2e9d284670ef676e177879990a4da0b568331deb))
* Update readme ([11d0942](https://github.com/3d-omics/mt_quant/commit/11d0942d89b37f5dc074fa37f8bcc2c2f7b9a29a))
* Update readme ([19c9962](https://github.com/3d-omics/mt_quant/commit/19c9962954ac761f07e70195353abd23893b8ddd))
* Update readme ([a6e5095](https://github.com/3d-omics/mt_quant/commit/a6e5095edb1fc0bdddf980097b7c64e5ba3bcfb1))


### Bug Fixes

* add prefixes to htseq and subread outputs to prevent a prefix loop ([974402a](https://github.com/3d-omics/mt_quant/commit/974402af461a49b32e796725ecc51129874caa8b))
* add user and permissions to rsync ([b52d0fd](https://github.com/3d-omics/mt_quant/commit/b52d0fd3a81ee0adb3d17830879e6b0cf5f2cf18))
* **folder:** consistent folder names ([89603a6](https://github.com/3d-omics/mt_quant/commit/89603a66b069a4a78d24cee29d251698153ff33a))
* **folder:** consistent folder names ([22f9066](https://github.com/3d-omics/mt_quant/commit/22f906638e00b7709d38ce9721f0e94faee23569))
* **folder:** consistent folder names ([15ba09b](https://github.com/3d-omics/mt_quant/commit/15ba09b84c9261807c1470123660efeda75889f8))
* **folder:** consistent folder names ([446d25b](https://github.com/3d-omics/mt_quant/commit/446d25bf6a1a97e9bf227285e05b42f6cddd16cc))
* **gha:** add temporary snakemake-action ([78945ad](https://github.com/3d-omics/mt_quant/commit/78945add1fe942dde7918529a1f8f7396a282c13))
* **kraken2:** use for loop over mapfile ([7a6afbf](https://github.com/3d-omics/mt_quant/commit/7a6afbfcada8a51542bc48189b206d868261e02b))
* **kraken2:** use for loop over mapfile ([25faa7b](https://github.com/3d-omics/mt_quant/commit/25faa7bbb4fa8f1f3354f1734c280f152a5ac81e))
* point to the correct branch in the badge ([aa9c2b6](https://github.com/3d-omics/mt_quant/commit/aa9c2b6bbb3bd67f2da2d7e380367cd6eff82934))
* **pre-commit:** update hooks ([b04c596](https://github.com/3d-omics/mt_quant/commit/b04c596b9a81aa46656ec165b0e133278d05e4c5))
* remove unused code ([8f1b473](https://github.com/3d-omics/mt_quant/commit/8f1b4734df885ed1f4a53a08eab8edced510181d))
* **report:** raise preprocess report ram ([80c84bc](https://github.com/3d-omics/mt_quant/commit/80c84bcd31235b7181c633e03933a67f04811f7c))
* **report:** raise preprocess report ram ([4a7f1b7](https://github.com/3d-omics/mt_quant/commit/4a7f1b73500aa95718df6832d8a9f992a385dd26))
* **report:** rename rules ([0a2197b](https://github.com/3d-omics/mt_quant/commit/0a2197bcfa49d2b45ee41f9c4df1d504c03cbbf5))


### Miscellaneous Chores

* release 1.0.0 ([0a18829](https://github.com/3d-omics/mt_quant/commit/0a1882919a09db7cdfb57f5d014bb15ab5c8dba3))
