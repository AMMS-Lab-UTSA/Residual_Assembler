# History rewrite of 2026-11-07

On 2026-11-07 the complete git history of this repository (Residual_Assembler)
and of its companion, [UMAT_source_transformation](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation)
(UMAT-OTI), was rewritten and force-pushed. **Every commit id changed.** Every
branch and tag was rewritten, so none of them carries an old id. Old objects can
stay reachable for a while through the pull-request references of pull request
#1, through forks and through caches, until GitHub removes them; their removal
is being arranged with GitHub Support. This page explains why, how it was done,
and how to translate an old id to its new one.

## Why

Both histories contained third-party code without a licence that allows
redistribution:

- the Fortran examples and input decks that accompany the textbook by Dunne and
  Petrinic, *Introduction to Computational Plasticity* (Oxford University Press),
  in both repositories;
- a single-crystal plasticity UMAT after Huang, as modified by Kysar, with its
  documentation, in this repository.

Removing the files from the current branches was not enough, because every old
commit still served them. The rewrite removes them from every commit.

## Method

- Tool: [git-filter-repo](https://github.com/newren/git-filter-repo) v2.47.0.
- Every branch and tag of both repositories was rewritten. (This repository has
  no tags.)
- The paths listed below were deleted from every commit, and the same file
  contents were also stripped by blob id wherever else they occurred.
- Empty or degenerate commits were **not** pruned, so the old-to-new map is one
  to one: each of the 154 old commits of this repository has exactly one new
  commit, with the same author, committer, dates and message. Abbreviated
  commit ids quoted in commit messages were translated to the new ids (16
  messages in this repository).
- Apart from the deleted paths, the trees of this repository are unchanged. The
  tree of UMAT-OTI tag `v1.2.1` is identical before and after, so its transform
  fingerprint is still `80a2015b780b1f25` (harness `56e3725136cb897e`).
- The rewrite is deterministic: re-running it on fresh copies gave the same new
  ids.

### Side effects

- **Signatures are lost.** A rewritten commit cannot keep a signature made over
  its old content. In this repository one commit was signed: the merge commit
  GitHub made for pull request #1 (old `009ca3a`, new `45a9212`). It is now
  unsigned. In UMAT-OTI one signed commit (old `d817003`) is likewise unsigned.
- **Line endings (UMAT-OTI only).** No file of this repository changed except by
  deletion. In UMAT-OTI, old versions of evidence files that quoted the textbook
  were replaced by redacted versions: where the project had already published a
  redaction, that version was used; nine older versions were redacted by
  applying the same text substitutions by rule; and for one frozen snapshot the
  manifest and checksum list were recomputed for its redacted file. The redacted
  versions of the two `convergence_rows.csv` files use
  LF line endings where the old versions used CRLF, so those old commits now
  show LF (CRLF to LF normalisation).

## Removed paths

### Residual_Assembler (this repository)

Deleted from every commit (39 files in all):

- `UMATs/ICP/elasticity/` (3 files)
- `UMATs/ICP/plasticity_exp/` (8 files)
- `UMATs/ICP/plasticity_imp/` (8 files)
- `UMATs/ICP/spin/` (6 files)
- `UMATs/ICP/visco/` (8 files)
- `sources/license-unknown/Huang_Kysar_Single_Crystal_UMAT/` (5 files)
- `UMATs/ICP/index.html`

The other files in `UMATs/ICP/` (the `UMAT_*.for` routines of the ABAQUS-US
collection, MIT licence) are not affected. Before the rewrite all 39 files were
still present at the tip of `main` and of most other branches; only the tips of
this release branch and of `cross-platform-hardening` had already dropped them.

### UMAT_source_transformation (companion)

Deleted from every commit (43 files in all):

- `UMATs/UMATs/ICP/elasticity/` (3 files)
- `UMATs/UMATs/ICP/plasticity_exp/` (8 files)
- `UMATs/UMATs/ICP/plasticity_imp/` (8 files)
- `UMATs/UMATs/ICP/spin/` (6 files)
- `UMATs/UMATs/ICP/visco/` (8 files)
- `paper_results/actual_umat_higher_order/code_imp/code_imp_oti_combined.f90`
- `paper_results/actual_umat_higher_order/code_imp/code_imp_oti.f`
- `paper_results/arc_791506/paired_batch/validation/code_exp/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/code_imp/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/elastic/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/spin_elas_def/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/spin_elastic/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/visco_beam/original_umat_with_ddsdde_sdv.f`
- `paper_results/arc_791506/paired_batch/validation/visco_imp/original_umat_with_ddsdde_sdv.f`
- `UMATs/UMATs/ICP/index.html`

## UMAT-OTI commits this repository relies on

| What | Old id | New id |
|---|---|---|
| UMAT-OTI 1.2.1 (tag `v1.2.1`): the `bridge` pin in `pyproject.toml` | `bf3f029d9711a2ac28fd011b6b373dd64ae24fb4` | `bb099d85b080a295a443773c823305671164653f` |
| UMAT-OTI audit baseline (`tools/audit_integrity.py`) | `3f340b72ec4d9d1f474fdb2e0f6bda1c41bfaa76` | `9f1f9d031f765fa81830ebf1140d7edd38d1a64e` |
| Residual_Assembler audit baseline (`tools/audit_integrity.py`) | `0789984c6c2750091d0f43b9ba42b84b3d890556` | `24d6c8fb743f4acf36b262e1aeecd12bf9599c83` |
| Previous `bridge` pin (UMAT-OTI) | `faa0360` | `95b27ff` |
| Residual_Assembler `main` at the time of the rewrite | `30a5ab4d0086` | `6e4c321362f4` |

## Commit ids cited in this repository

Commit ids cited in this repository's documents, evidence files, code comments
and tests were changed to the new ids, keeping each citation's abbreviation
length. Three kinds of ids were deliberately left as they were:

- the directory names of the frozen UMAT-OTI evidence snapshots, such as
  `umat-064583d71b90_ra-f928acec970a`: they name directories that still exist
  under those names, so they are labels, not citations;
- ids of local working commits that were never published (for example in the
  second column of the tables in `docs/BRANCH_IMPLEMENTATION_AUDIT.md`): they
  were never part of either public history, so the rewrite did not touch them,
  and they could not be resolved on GitHub before or after;
- captured and proof-hashed evidence records keep their pre-rewrite ids, which
  are looked up through the map below: the raw clean-clone report and its
  milestones (`docs/evidence/recovery_clean_clone.json`), the proof blocks and
  captured gate heads of the usage and portability runs
  (`docs/evidence/usage_examples.json`, `docs/evidence/recovery_portability.json`),
  and the prose and captured output that quote the same runs
  (`docs/evidence/recovery_clean_clone.md`, `docs/evidence/recovery_integration.md`).
  Their sha256 digests describe the bytes as captured, so the ids stay as captured.

The pre-rewrite ids kept in those records:

| Repository | Old id (kept) | New id |
|---|---|---|
| Residual_Assembler | `5ea3da4bdd250fe6eea79677cec44e59e86de197` | `54876c47e2ce22b2562a20d7c813629dabebd071` |
| Residual_Assembler | `69ac77cf4865ef5cfc2a46fad987bb7fb6d73e21` | `5cfa1c85fe8a8e7cee8c6eb89c7577775ee22d68` |
| UMAT-OTI | `dbd8ff229023b57e4231a2250a54e97dc111f547` | `397f53a47e6f2aa6a517c7eebd46254cdf094c41` |
| UMAT-OTI | `d34c7fff51e13285efaf13a57917fa34034f57cb` | `7eb31b734e1c8d8601770a20981a29634dd9f7e8` |

## Translating an old id

Look the old id, or its first seven characters, up in the table below. The
UMAT-OTI ids this repository relies on or keeps are in the tables above.

## Residual_Assembler commit map

All 154 commits, oldest first (author date).

| Date | Old id | New id |
|---|---|---|
| 2026-07-13 | `c01202978f8a81b308b5c26d260359507ed21933` | `57be11683b5bd8e1fa2e0663ae9ae8133e9cd0d8` |
| 2026-07-13 | `00c784adbe90eec4d5ac0a6349f39742b5a47176` | `69e3a68b36fc8442c1d90ac8687893ec52dd5daa` |
| 2026-07-27 | `4a92c508f6b2fba517db95abbcb94e624cd6e79b` | `7f5f63b8f56bac09c18b06d05ad2e2a15baa684e` |
| 2026-07-31 | `bfde4d027652b1b60d306a42bdc6746f64ddb24f` | `ed6103c49a0a662e826a5bde96abbb576171247e` |
| 2026-07-31 | `009ca3aad4429cc09f7623926641228883bde614` | `45a9212b6ea7d1ac193a484067f047e22089f86b` |
| 2026-08-21 | `49393a1a7cd96d41ab833b21b386f5aaf25dc56a` | `7c68b3e7c301b1d97d04276bcbf242987384f72d` |
| 2026-08-21 | `4739b1b8cb14793c2c1bf460ad2e806c451112d8` | `0fb19c15c3af99d58428f7ab0b1f96babfa72c85` |
| 2026-08-21 | `838d66bc23cd23b3287611f4dde8f57d6bfe3902` | `5efdd0cd7955b376cb24c9f4d64fd0fcdc603b26` |
| 2026-08-24 | `949e0dc3a0af343be1e51b33ece1bda253064bcd` | `09c69a6b54f3d3f7cbf835c4898141c3620375c9` |
| 2026-08-24 | `15b4d6183f2aa55eca42fb32a54e14e347ae243b` | `41bca8a45a25e236b9b6c03412eee4cc0e445601` |
| 2026-08-24 | `74deacc28db176696d50deee78e16eee81b19402` | `86f497ea005e5462920546a11a72e3cdc0c86d0a` |
| 2026-08-25 | `4b76c3cc6100873904cff688e3b81ecfd39339b2` | `5795e2c0ed47f0d4809773c6717c6bb385f2fd90` |
| 2026-08-25 | `9694033c50efe630de26b4a253ea58295be6c29a` | `783ea721f0e1527162fc84cc56ed1c69643cdfce` |
| 2026-08-25 | `cec69908e79e8d87a05af6fae8e391117bda5bba` | `92288df9538d3c8cc62fa11013c970180b96115e` |
| 2026-08-25 | `3b7411e4ede328063b38b2537eb3db8612e5bdda` | `4fb873a747ca2b486e9dc535a136112dc32ca31a` |
| 2026-08-25 | `00b715663e41661735987ab92e02b117261f6a1c` | `3634be5045e4ea9ea7ab036c88758877fb2c787b` |
| 2026-08-25 | `80fd4afddd1cdd4c2807baddc44cc616dfe7f09d` | `2b0496e2eb9898b1d2d7dc03764494ddc6e9a3ad` |
| 2026-08-25 | `3802d93012620efef5cfd3a9376fd44d11e03893` | `3680e1a7e0bee2bb123ea186c12d5b4b3dc2f727` |
| 2026-08-25 | `30dcdfc02c8fe1857ed189bdd450cab2a2007433` | `3b1f1dd851ab7814b95d946e77c79727b6c28eb7` |
| 2026-08-25 | `f928acec970acf68670eba9102e6b962dabb094c` | `7176650dcf693bae548457d428ac094d3e06c505` |
| 2026-09-10 | `c7d4fb6d30d36767c60bd41c49914f05bbebd104` | `401bc76c5ba41193331f9cb1f416a39f54c8a3a0` |
| 2026-09-10 | `031fcb75508e25bab52ea377172e1bcdbb3d86ec` | `f2f5307af62690c399634eaf495453b6b090e25c` |
| 2026-09-11 | `c03c24c1d686ca2d2b21e39e674615f40472f137` | `d7a3e715f01e5e3bfc6471a4aacd0a08d8e5ab18` |
| 2026-09-11 | `9e8ecb51780e28ea85480a23ecea85e34bc35c0f` | `989e63db1a432946dc3726f311d9770568c62e8b` |
| 2026-09-13 | `bfd4302c047cfab2c6931203ce8f157d1b5bb2c2` | `f04b85f529088db99c840ff0c9867c58050683cb` |
| 2026-09-14 | `ae3b804f094e968dea046398a43d749f57a277d9` | `7a0f55f6ed7b94530fb520685fb6da71c44acec1` |
| 2026-09-14 | `75b7279e4467847469eb11f4c21e532397cad093` | `622e12d22d9f1cb51e50a6c109f65bab3ec68a39` |
| 2026-09-14 | `8902bef279d3758a31bda5b59e4a725f35bda037` | `a9ac801dfb3e0cb6154d3e80fbbc9d39b9e5be18` |
| 2026-09-14 | `b9e4ee5477cc9a72eadb1f5d9a4053c5f1fd4522` | `6eef5c49765066e5ed3526363026129a5be6dc04` |
| 2026-09-14 | `475729bc87f4f1acd08cd231070897b5b9d80011` | `832a150bd246260e59ced61820a84d666c024374` |
| 2026-09-14 | `14abe1d5cd7fd6915eea24e8ee0b23b22d35571b` | `cef05257ed7f060cd390f6ab648900ef50e90b18` |
| 2026-09-14 | `4c9c25747ea0469f6fd9fa58e57940f412fa7e53` | `6cad14ecd81d24133d30141c301ee6d9db452e81` |
| 2026-09-14 | `e45a6a85571da48313ca63c1a6fce3024422ac80` | `4f42c01fe387ac1f3e91b4e7f1406f48e10f5fcb` |
| 2026-09-14 | `004abc117d6ff47fd4e2fb9ebad4e65d00116546` | `d6d73b43747a4d1cb3d861583996724415232cc9` |
| 2026-09-15 | `ea6c7ef33cf2dc33a930f51938de1dd148721b13` | `1823dfe91f13371807674a3e539675e856bacc52` |
| 2026-09-15 | `0435eab5c092821ad9c8f5a548f7409255af5df1` | `2c809c9d8720ece1985ca65e4aab19972b2792cc` |
| 2026-09-16 | `ab0ee10c6b842bb80fba189279a9afb572766252` | `1748d45865add7a21a46c8e686b0aec9a46bda13` |
| 2026-09-17 | `053efc61abf1c975820eeff37aaedc9dae55e8fb` | `b46c44d37424c4b355b7c5531652e1aaf058f6bd` |
| 2026-09-18 | `f1d9a32caf43f6370e709959d7c51cd62097a0ee` | `e4c7b61e144b85e1bf96aebb3ea43f4ff429c3fc` |
| 2026-09-18 | `69ac77cf4865ef5cfc2a46fad987bb7fb6d73e21` | `5cfa1c85fe8a8e7cee8c6eb89c7577775ee22d68` |
| 2026-09-18 | `5ea3da4bdd250fe6eea79677cec44e59e86de197` | `54876c47e2ce22b2562a20d7c813629dabebd071` |
| 2026-09-18 | `0789984c6c2750091d0f43b9ba42b84b3d890556` | `24d6c8fb743f4acf36b262e1aeecd12bf9599c83` |
| 2026-09-18 | `0a680d92255ae0a0fa17ef5538dc19e6c97ec30e` | `7addeef86aad9158d0de674ed672ded35c77e5a8` |
| 2026-09-18 | `9dd6a8c04f6c610156e2716fd4d4934b65f2a07e` | `96856cfce8ed4c7aaf977d412c619ae2c90251df` |
| 2026-09-18 | `a7271cde2ede05da189df3ac8be4fe671db86425` | `d902c212f1594b70fca5caa78d966c8ca1daf854` |
| 2026-09-18 | `107e8c0fb7b19a4a6ecc1fca2919e99e1d0df248` | `2deb051f8da8b4b7f03a3bcff6f07385b3b3bc9f` |
| 2026-09-18 | `15a524014b89cff1c8ec35d03b00bc79cba03314` | `0beed34fe9301106bedea74c9043b06e2c9fd3b0` |
| 2026-09-18 | `b891d74fe6283eb7edc9f708934eeaadf0b5fc54` | `18f837e540267ab351862d921a3285c4a475fa51` |
| 2026-09-18 | `0b3746aa56d32d7dc8311f1fc0b2562437cd2984` | `f1542f08f13684fbc7612c05ad6e8b7a4c90b10d` |
| 2026-09-18 | `58cfcfdba4f1f75c1b6b252d300d3d958899f7cb` | `6ce545db2cae68a3be878b301c0f5787a0b5d54d` |
| 2026-09-18 | `7b396d6a148ad98ade2f0d6fd443226dfab8dd1d` | `ea71f639a1cebf5c843c10d28ffbf6a4cbd4db60` |
| 2026-09-18 | `be00f293c66433b79034baa98f7461c9aa9e845d` | `551901e9af722d6117d8d23ce9cb789ca7cf1143` |
| 2026-09-18 | `cb03bbb123e0ecd6546c5c7163897fdc59710e91` | `8d47bd9f56d7a77c1930987d971beea3ac8191c1` |
| 2026-09-18 | `c51b2a1be674d22b79e8c656d56b443ceed64893` | `59c1c038f2ec2cf36983e7c3593300f84bf3ef46` |
| 2026-09-18 | `e15beb9ef61a5a508234d06e09d5b105e83722e1` | `9da43eced440bb324c07f7b050fbb1eeeed4eca1` |
| 2026-09-18 | `e0003452636ace9b708f7db514bacd8318c2f5a2` | `572d036b2d8d985743c5f3ecf6e31ba283e9f3b8` |
| 2026-09-18 | `2210eebdc83d5871333344ffd812cf04d965db14` | `0b17e5788d8430fea9c42a5284aa2a5f17448e58` |
| 2026-09-18 | `b7a257d2486f24cfef9d4db36bcb785395e4d669` | `e34907ddcb2cd8602a967caff10ca373a8deb66a` |
| 2026-09-18 | `718ac7c404b5a296092223b20537d2d388979f49` | `b859cadf4adfe793ca95d10edf0948827ef7bf18` |
| 2026-09-18 | `5e49cc96b70684e330c38e1a4c5d298a7153bb40` | `d31511eea98d0a7dc348528b4f2824a398feb322` |
| 2026-09-18 | `9999d1229e05ba0ddb262b30470d35ba1f0ecadb` | `078144d45ec24abed799bd2f296c2e44f40124e9` |
| 2026-09-18 | `f2abb0b2f65e0b98bd4fbff4d59e519e4a5142d5` | `11f502ddcf052d342d4abdada24b5c7a2881bcf9` |
| 2026-09-18 | `e195c085b677a7ae47c8a18f9019d054b2bfa7f0` | `7d35b29658579d8718e2cdf91ce764ad41fb5505` |
| 2026-09-18 | `3504a0294107e17d32de451f8dd3b7a94abccd1c` | `edfb94c268428e46ea9a792ca91cbc93c1ec0da0` |
| 2026-09-18 | `e03ef6054b6cc4201563d207bf6c6610d9289095` | `cd3e815a236787e150b1d1dad5fbc4fc5df21903` |
| 2026-09-18 | `8daf75fbfed6edb2ce46abddad19ae5d95408efd` | `202bacb6df55fe23dae624eeed867b4ab1956576` |
| 2026-09-18 | `71cfe8660d1eaca87fd08fb895916efe9ae0ca03` | `e7d3cba9903cd0774cc6e24b21ef762518bf1440` |
| 2026-09-18 | `1bc2009a1f45584eb4693167acab7ffa5c90bb8e` | `84a20ad145faf9074ac733d433f0a182d4954f9d` |
| 2026-09-18 | `924ee8cda189bd517170687f237aa4e9b4313ccc` | `d05d6b043fc4ea1dcfaf1ec718d3e5324889fee9` |
| 2026-09-18 | `ea1ed80f06a5bee78b1ecb1572255de7bcca86d7` | `c98b10718d62e46c925fc0579146f470418a0a70` |
| 2026-09-18 | `133d8c897f9f73666de7bc968e7b2b538fec3a23` | `f120f81a578109ce576d0698b4584111266e5750` |
| 2026-09-18 | `d27a6a0fc428440a0c83283faaa23b327be25714` | `8b4c805e74afbd2655792a15d1308947880a67a8` |
| 2026-09-18 | `3bab862fe10722522521ee8b5fa952f38846a1b4` | `b3442e6362b1b786ad00144f035b66aa2b5a5669` |
| 2026-09-18 | `a3916bec1547775ccc489195367df8799553c4e3` | `fa52464dbed92a27dd7923b034c1d55998ffc23b` |
| 2026-09-18 | `77743136353773d40cb0fbd522914bcfbb3eaec7` | `5b7866793d2565f785427041ac475077c48d0ae1` |
| 2026-09-18 | `de8ec64a34a61e442abb2c653a2c5a561b5dffd4` | `085cd999f0a1bc352af0b48372b8477242d3dd3e` |
| 2026-09-18 | `4d34c746ea255286f02931d4944a6082c30fd98c` | `992c517db0ab790e06ed66fe7ac1df60de010048` |
| 2026-09-18 | `9e5e3b9f0cbbcf02229443fca147f036cc944d6f` | `d42498d4470ff494c602e220c6f63b6eb7075221` |
| 2026-09-18 | `de928908cac269e94e0ae8c98d3ff8333f74be00` | `0267ff530d5d21b9e499514aeebb2a906f0ed3c0` |
| 2026-09-18 | `07703d3345a7487c33e33ac61f238420875c82ae` | `e77695755a9ba9ce9c1114a1a97475e078f73bcc` |
| 2026-09-18 | `6df72c189390cca7d90df5d777ccf2b56175e9dc` | `b742a05646ff4a571938e574ba38464e6cb61986` |
| 2026-09-18 | `fb3d059b8d73ae14e4c5db52906256581d1b6fb2` | `28576b0ac8bfdabfd9e1f63347c33af3a1a70318` |
| 2026-09-18 | `a392cd5b1277cce8c03093f7ac51ee1754dddf3a` | `ea2baee5f4d18a94b8521fcf82cdf7747a2120f8` |
| 2026-09-18 | `d211658f29e77cb35bce5526248974e5b01b9413` | `ecff15ba5b37ff334a84e5b9fee36189dd71a0dc` |
| 2026-09-19 | `7dfa511fa1c06497cc29af236b5c037b129cf4db` | `22dd957ee3e354576db93f9da203cf452985ee68` |
| 2026-09-19 | `8e4cf2ff660ebeda86ffab13359ccf0479895b8f` | `a19127b80f49b9b9deee22d95a80a2c3d806c9f5` |
| 2026-09-19 | `d3cdfd51503fd30e4c1c575acdade43a2fecd9d0` | `4d7bbf7817e70ee1b325bc22461b919b79d2ebe9` |
| 2026-09-19 | `60c8d44a2fafa936e0d855e47a8e57860a1e7a40` | `f013577554e851e94e844f47b58202fcd0baf1d3` |
| 2026-09-19 | `a288783a0965906eee9e9e123117af5fa8cc9511` | `c9bb66373d5f05cc0a3a7238dfbf324d937dbb3d` |
| 2026-09-19 | `83fc0a604c32e59ee67267c4463a2f1ea675f539` | `1ac4d68528ee63d6ecb6e838d3ec2c23192a5b16` |
| 2026-09-19 | `7099a1a8d0937de6f74dec1ab552b852ae34c967` | `b34a96e1270d712233cf182dd96b249655a430a1` |
| 2026-09-19 | `fe85ab84748b2c03d20d91d30e082b6f2d3cab75` | `ca043008ba5b5c639cdb6646259f3a82480db505` |
| 2026-09-19 | `225ad6197609e889247db408ee36d7ae6326a43a` | `07582aed7a4b9c3aec487d0e7de8beada0b2e4df` |
| 2026-09-19 | `f4a2073e916d3d1595dc366fffe0e5a3cbea9127` | `83466d6e65cdaf229845271ade268bf5817ca85b` |
| 2026-09-19 | `b8f7a85e4510fe74221cf91b31b70b1997063f54` | `c519df42ad1b3f38f932c03475ddd1002937170c` |
| 2026-09-19 | `f2f326db41b533916a123170a2015de415f05a0a` | `5995d3f06fb69bb973bacbc27400f9328711902b` |
| 2026-09-19 | `bf3600cb58d8f267020573476dc26c826c5048f5` | `db64d746fb0176ad94100574ac94efa64909af30` |
| 2026-09-19 | `acf7c207f788a43fcc9ac62109fcdc1c2cb224a4` | `42336a1696c5d5abc5ad64396c56f80c3c1ffbe6` |
| 2026-09-19 | `83168037cefc89487046915bd8f8177932449725` | `8ab62951008f04228c5c84e82a4b2d4de54e3917` |
| 2026-09-19 | `5694af029d2269f4eb6579e25af727c55efbc29f` | `64cbbac8ac79a1e34b3c56369bd86142a41d2122` |
| 2026-09-19 | `24a0f13759afa77d33d2b657f2d69e9b6f2ecd2f` | `c9bbeb8e0e009678285ac13197ff29d1d67c7804` |
| 2026-09-25 | `30a5ab4d00867bc419ae8e8e866530fd7e2e9e65` | `6e4c321362f464c51f9dbbbba510ff5bc618257d` |
| 2026-09-30 | `d74bce815a7f52b93bf982496a4bb744f7bdd829` | `39c10b337a87b53e048aff6e93530da61860da19` |
| 2026-10-01 | `0c48c8891300081ce6bc93014eae58393d272edc` | `e39b505d113cc3177ce24a4f774e2dc2bd63db0c` |
| 2026-10-02 | `c8d0b0ab10bda731fdfd3fc14c9740e7264c03b3` | `899977fdf587c9f22b520f130cdd051368f539e9` |
| 2026-10-02 | `7fcfa3f853333dc9aa6d74927262c49d5c2c83a9` | `2c9323857ab73ec6161acce3ed08c6c3d5ee6421` |
| 2026-10-02 | `6b5357bbf2ef650b2709dc5cb860e1ed9cbe448b` | `f62a6bfe1f7ed916c8404abf8c5a0e42005e9137` |
| 2026-10-02 | `681f71c57df0bfa6c8814a6f553aa67b0b3ab3b9` | `ef0863f95e98e5b23fdab08a81572792dd88555a` |
| 2026-10-02 | `7a490d3220868ea40bc57ba23e50f406df765423` | `6a5219358745a132eb8fa91b529c78baa78a6c2d` |
| 2026-10-02 | `c0dde1399dd188cfa31ceb1312f6bac743af408d` | `d68aece19be5d03de1c20dff222ee464628d3305` |
| 2026-10-02 | `82a7d3b1b9efe7b1bfab2e5964e001ac77083dae` | `2fe044956cc2aad24bd331abfe2031f38d681005` |
| 2026-10-02 | `644041e7c17ef503ff9658d21de81783f18983a9` | `df557d5b5710730713413f2734107c6a36d7dc24` |
| 2026-10-02 | `cde44a0a80ce3e9bc534ca9e8df748f0b3ec60e0` | `6066af8a1d130baec339dc9a9ba00f51e5b55213` |
| 2026-10-02 | `2712a39d5bfb9086bdd524f562ba60c683521138` | `63a9896754a4fe1ab46a77c6a388a4916431a0c8` |
| 2026-10-02 | `ffd986ef437a158652e43411b28ec51dd6479aa0` | `39554ad04adb00350f3c5ffb5c6f840aedfee532` |
| 2026-10-02 | `5241628418f8bc7702ccf78a9ef7dcb3083b4deb` | `c2d85fe99003cbc96d4a75c702b9b66d0c2c9a2e` |
| 2026-10-03 | `748ca605999add2f343d27e58ba6b8c0cf00fbce` | `e90a5505ff38907e1da4dead4848ea61e18dd511` |
| 2026-10-03 | `a7ccce6604d57f25e485614ca69af4a02aff6c59` | `a86ca58be7fffe01d7c1127008dd8ae2d4d7908f` |
| 2026-10-03 | `a5faf5ce334e223ceb9eeffdb5e06042f34aeb2a` | `ebcae0149765dad01f6ace235b49d367fdd1be6c` |
| 2026-10-03 | `b5243451b2d8add9766834db60a4516a543081f4` | `4a9c9c6aa6b8373825ac2e57fcd4bb13fe1b9301` |
| 2026-10-03 | `f1abeb991bb9de3710977d9a5359b7e257f030bd` | `5376416b2d3c058ac7bcfe54fdf3ca3c7240f6af` |
| 2026-10-03 | `a6b0b4180cb45390c77a003f139ffcb394b1e201` | `21b16c9ce97ee3b932f817ee81fff93b63a1ced8` |
| 2026-10-03 | `74eb0755ff066044c05c2889e50aaa43b0873da2` | `cab00df118148e980dc8d9a29ebaf4c9de228c8d` |
| 2026-10-03 | `1527a0c6baeb176680330789bb09744e58f1a64e` | `7edf67e29728b7bda512683a293482fab6b5a1ce` |
| 2026-10-03 | `df6a86f28a72af48c606b9a51e7a02fee44f5ada` | `147764c3b51195365f8942512b5d4cb9bac73ac5` |
| 2026-10-03 | `808fcb6aa46c560d5c89db0172875f168ec920a1` | `814edb272182c9d2d36489e038532b11b821beb5` |
| 2026-10-03 | `bcda33bd109983bd5812ee2b4f0f913dbe9281ef` | `fb15c637a7e99b6b2a302d55cb0be2430ab9a890` |
| 2026-10-06 | `8b9613b157742b7589396d8c18ee6dc2319130bd` | `94675328534c88d5a09692c80dfbda24f4f6f359` |
| 2026-10-06 | `a353eead760895957a39e38e193c2817f3ad4c78` | `314c2204c697eb4356bf618615aedf9aff2b4b2c` |
| 2026-10-06 | `6441924f64aa12a2a16a8885e62081095075dd77` | `d259c209f6254886e997378e1ad8bdfedc4efb52` |
| 2026-10-06 | `1b663ba0577cb2e4fac74a548748c07f97bc3b6a` | `13a6904ae095fda4285db276d182c9d06e142d24` |
| 2026-10-06 | `811b13c3e4e04af972246db5e460cfcecca1ba6d` | `4cb0969b8fd82fe73709337011da8506ed57b3b1` |
| 2026-10-06 | `42c5e0d1a11844cfdf521b6fdc0d66fd42adc580` | `eb7cf0403750456f9d17d2f2d15815ba55fcb995` |
| 2026-10-06 | `c88a1dee96cf7e3956c316a2a23dd568c2805b81` | `5abb2b681b9cf6066566b8e274e5cfa8d4065372` |
| 2026-10-06 | `2df8d5b10801c801048e0fc5fa12ac50bdc13c20` | `f38f5547ab84e58d4b5aea67148afd0ea6ba8004` |
| 2026-10-06 | `a64aff5eee5775ae731945ad949ba6c7db63d63c` | `51702240dc102ff114183ffe6bba333dc69d0f15` |
| 2026-10-06 | `a4e31e5aeb9aa653dade6cc3749ceacc0e6620dd` | `cb930d80d903e181927d9b53ccb92c5b6a658e1b` |
| 2026-10-06 | `f07bc29fe178b7171a14e0b5642e0907690ecd6f` | `f3efc659387139b84b8ff27f91ff676368856333` |
| 2026-10-06 | `582315097c0c0f267b29dd492b923e32133091f9` | `1bc377edc87a1fa91f5b7a7db3ed112274c8ddf5` |
| 2026-10-06 | `e906569fd8828c16b16de1ce6f3f26bc8618a3ca` | `b6473b7fc879046406f23380cf6fd152583cd6bc` |
| 2026-10-06 | `769db2dba0e4ca59a20668136fb7c521809d04e1` | `d6583bf0119c42b172e46a179235e5229893f432` |
| 2026-10-08 | `87205e31d027817639c504fd22611034c4fa0a98` | `e6523d61c8056cbcc1353635fcf2fcf5a6de8206` |
| 2026-10-08 | `65977a9576a55a2d3ff35d7166b115e79807b28e` | `cb45afc6706f305ea69a6555a9b75b9cbfd302e4` |
| 2026-10-08 | `412b68ef9b64b5e8a0fc683470f1f08a2bea72a0` | `07ce010e72d74c62a35d6c796e72a20ec4140738` |
| 2026-10-08 | `de7584e600fa4eb3852b2c787d48bb0b07ae5cb2` | `6a0d969770f367b55a5c0427af9d4f45c517c61a` |
| 2026-10-15 | `0f1e5ad4971e98eae59c52a70258c1586b563967` | `8466cf475c98a14a8d5d6271f467ba74077a0721` |
| 2026-10-15 | `e47090e0bbe0babe3e2bcbc6cb10f3d2095e8dfc` | `aabfa59f19c7bf9c9f06a4d2b047ff9cef9953b2` |
| 2026-10-15 | `0f9fcd9973a1a56062bc77ade105337c803a0fb3` | `0d8b10880073ec127a8172a8cff10763daceb97e` |
| 2026-11-06 | `3895fa4f066f03b7267c0801516935cf04909d32` | `b7cf10226416ba601076c41603eecf89ddc87029` |
| 2026-11-06 | `efa0418876ea33f21555d388b75bfefd4e5fe06f` | `446d0d73f9cf1ec09bf369480756f35c8804f9d2` |
| 2026-11-06 | `daa636ebd552567ada51306a9fa9d11de9edba67` | `b961853b442e797c9f252d3dd31bfddbc0bdf200` |
| 2026-11-06 | `44f6f90ad7454f51aa480a6518ed36212a09c88a` | `a004ad67f72af30d0da2c8850f2d9ee22a683d3b` |
| 2026-11-06 | `de86038d7dac58c655bc38d425cb764804e0a0d4` | `2341073fef96884f8a36ade936321b2a40c789ff` |
| 2026-11-06 | `a8056ca73ccd67b7a609444389f437cdc07e7001` | `0d94912ec7018f074debcae0c9dfd3d5836f2d50` |
