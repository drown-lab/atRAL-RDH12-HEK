$ErrorActionPreference = "Stop"

$outputDir = "Proteomic_Figs/Supplemental_Acute_RDH12_GO_Reactome"
$outputSvg = Join-Path $outputDir "acute_RDH12_5hr_GO_Reactome_supplemental.svg"
$maxTermsPerPanel = 16

$panelSpecs = @(
  @{
    Panel = "A"
    File = "GOanalysis/output/clusterProfiler_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_down/enrichGO_BP_manual_filtered.csv"
  },
  @{
    Panel = "B"
    File = "GOanalysis/output/clusterProfiler_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_down/enrichGO_BP_manual_filtered.csv"
  },
  @{
    Panel = "C"
    File = "PathwayAnalysis/output/Reactome_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_up/Reactome_manual_filtered.csv"
  },
  @{
    Panel = "D"
    File = "PathwayAnalysis/output/Reactome_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_up/Reactome_manual_filtered.csv"
  },
  @{
    Panel = "E"
    File = "PathwayAnalysis/output/Reactome_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_down/Reactome_manual_filtered.csv"
  },
  @{
    Panel = "F"
    File = "PathwayAnalysis/output/Reactome_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_down/Reactome_manual_filtered.csv"
  }
)

function Escape-Xml([string]$value) {
  if ($null -eq $value) {
    return ""
  }

  return [System.Security.SecurityElement]::Escape($value)
}

function Wrap-Text([string]$text, [int]$maxChars) {
  $words = $text -split "\s+"
  $lines = New-Object System.Collections.Generic.List[string]
  $line = ""

  foreach ($word in $words) {
    if ($line.Length -eq 0) {
      $line = $word
    } elseif (($line.Length + 1 + $word.Length) -le $maxChars) {
      $line = "$line $word"
    } else {
      $lines.Add($line)
      $line = $word
    }
  }

  if ($line.Length -gt 0) {
    $lines.Add($line)
  }

  if ($lines.Count -gt 2) {
    $secondLine = $lines[1]
    if ($secondLine.Length -gt ($maxChars - 3)) {
      $secondLine = $secondLine.Substring(0, $maxChars - 3)
    }

    return @($lines[0], "$secondLine...")
  }

  return $lines.ToArray()
}

function Interpolate-Plasma([double]$value, [double]$minValue, [double]$maxValue) {
  if ($maxValue -le $minValue) {
    $t = 0
  } else {
    $t = ($value - $minValue) / ($maxValue - $minValue)
  }

  $t = [Math]::Max(0, [Math]::Min(1, $t))

  # viridis plasma-like palette, reversed so lower adjusted p-values are yellow.
  $stops = @(
    @(240, 249, 33),
    @(252, 166, 54),
    @(225, 100, 98),
    @(177, 42, 144),
    @(106, 0, 168),
    @(13, 8, 135)
  )

  $position = $t * ($stops.Count - 1)
  $index = [Math]::Min([int][Math]::Floor($position), $stops.Count - 2)
  $local = $position - $index
  $red = [int]($stops[$index][0] + ($stops[$index + 1][0] - $stops[$index][0]) * $local)
  $green = [int]($stops[$index][1] + ($stops[$index + 1][1] - $stops[$index][1]) * $local)
  $blue = [int]($stops[$index][2] + ($stops[$index + 1][2] - $stops[$index][2]) * $local)

  return ("#{0:X2}{1:X2}{2:X2}" -f $red, $green, $blue)
}

function Get-Nice-Max([double]$value) {
  if ($value -le 5) {
    return [Math]::Ceiling($value)
  }

  return [Math]::Ceiling($value / 2) * 2
}

New-Item -ItemType Directory -Path $outputDir -Force | Out-Null

$panels = @()
foreach ($spec in $panelSpecs) {
  if (!(Test-Path $spec.File)) {
    throw "Missing input file: $($spec.File)"
  }

  $terms = Import-Csv -Path $spec.File |
    ForEach-Object {
      $foldValue = if ($_."FoldEnrichment_num") { $_."FoldEnrichment_num" } else { $_.FoldEnrichment }
      [pscustomobject]@{
        Term = $_.Description
        Fold = [double]$foldValue
        Fdr = [double]$_."p.adjust"
        Count = [double]$_.Count
      }
    } |
    Where-Object { $_.Term -and $_.Fold -and $_.Fdr -and $_.Count } |
    Sort-Object @{ Expression = "Fdr"; Ascending = $true }, @{ Expression = "Fold"; Descending = $true } |
    Select-Object -First $maxTermsPerPanel

  if ($terms.Count -eq 0) {
    throw "No plottable rows in $($spec.File)"
  }

  $panels += [pscustomobject]@{
    Spec = $spec
    Terms = @($terms)
  }
}

$svgWidth = 1420
$svgHeight = 1920
$panelWidth = 700
$panelHeight = 620
$leftMargin = 360
$plotWidth = 220
$plotHeight = 450

$svg = New-Object System.Text.StringBuilder
[void]$svg.AppendLine("<svg xmlns='http://www.w3.org/2000/svg' width='$svgWidth' height='$svgHeight' viewBox='0 0 $svgWidth $svgHeight'>")
[void]$svg.AppendLine("<defs>")
[void]$svg.AppendLine("<linearGradient id='plasma_adjust_gradient' x1='0%' y1='0%' x2='0%' y2='100%'>")
[void]$svg.AppendLine("<stop offset='0%' stop-color='#F0F921'/>")
[void]$svg.AppendLine("<stop offset='20%' stop-color='#FCA636'/>")
[void]$svg.AppendLine("<stop offset='40%' stop-color='#E16462'/>")
[void]$svg.AppendLine("<stop offset='60%' stop-color='#B12A90'/>")
[void]$svg.AppendLine("<stop offset='80%' stop-color='#6A00A8'/>")
[void]$svg.AppendLine("<stop offset='100%' stop-color='#0D0887'/>")
[void]$svg.AppendLine("</linearGradient>")
[void]$svg.AppendLine("</defs>")
[void]$svg.AppendLine("<rect width='100%' height='100%' fill='white'/>")
[void]$svg.AppendLine("<style>")
[void]$svg.AppendLine("text{font-family:Arial,Helvetica,sans-serif;fill:#000}")
[void]$svg.AppendLine(".panelLetter{font-weight:700;font-size:28px}")
[void]$svg.AppendLine(".axis{font-size:13px}")
[void]$svg.AppendLine(".term{font-size:11.5px}")
[void]$svg.AppendLine(".legend{font-size:12px}")
[void]$svg.AppendLine("</style>")

for ($panelIndex = 0; $panelIndex -lt $panels.Count; $panelIndex++) {
  $column = $panelIndex % 2
  $row = [Math]::Floor($panelIndex / 2)
  $originX = 24 + ($column * $panelWidth)
  $originY = 24 + ($row * $panelHeight)
  $plotX = $originX + $leftMargin
  $plotY = $originY + 66
  $axisY = $plotY + $plotHeight + 6
  $panel = $panels[$panelIndex]
  $terms = @($panel.Terms)
  $minFdr = ($terms | Measure-Object -Property Fdr -Minimum).Minimum
  $maxFdr = ($terms | Measure-Object -Property Fdr -Maximum).Maximum
  $minCount = ($terms | Measure-Object -Property Count -Minimum).Minimum
  $maxCount = ($terms | Measure-Object -Property Count -Maximum).Maximum
  $maxFold = ($terms | Measure-Object -Property Fold -Maximum).Maximum
  $xMax = Get-Nice-Max($maxFold * 1.08)

  [void]$svg.AppendLine("<g id='panel_$($panel.Spec.Panel)'>")
  [void]$svg.AppendLine("<text class='panelLetter' x='$originX' y='$($originY + 18)'>$($panel.Spec.Panel)</text>")
  [void]$svg.AppendLine("<line x1='$plotX' y1='$axisY' x2='$($plotX + $plotWidth)' y2='$axisY' stroke='#000' stroke-width='1'/>")
  [void]$svg.AppendLine("<line x1='$plotX' y1='$plotY' x2='$plotX' y2='$axisY' stroke='#000' stroke-width='1'/>")

  foreach ($tick in @(0, [Math]::Round($xMax / 2, 1), $xMax)) {
    $tickX = $plotX + ($tick / $xMax) * $plotWidth
    [void]$svg.AppendLine("<line x1='$tickX' y1='$axisY' x2='$tickX' y2='$($axisY + 4)' stroke='#000'/>")
    [void]$svg.AppendLine("<text class='axis' x='$tickX' y='$($axisY + 20)' text-anchor='middle'>$tick</text>")
  }

  [void]$svg.AppendLine("<text class='axis' x='$($plotX + ($plotWidth / 2))' y='$($axisY + 42)' text-anchor='middle'>Fold enrichment</text>")

  for ($termIndex = 0; $termIndex -lt $terms.Count; $termIndex++) {
    $term = $terms[$terms.Count - 1 - $termIndex]
    $y = $plotY + 22 + (($termIndex + 0.5) / $terms.Count) * ($plotHeight - 44)
    $x = $plotX + ($term.Fold / $xMax) * $plotWidth
    $countScale = if ($maxCount -le $minCount) { 0.5 } else { ($term.Count - $minCount) / ($maxCount - $minCount) }
    $radius = 3 + (9 * [Math]::Sqrt([Math]::Max(0, $countScale)))
    $color = Interpolate-Plasma $term.Fdr $minFdr $maxFdr
    $wrappedLines = @(Wrap-Text $term.Term 48)
    $labelX = $plotX - 10
    $labelY = $y - (($wrappedLines.Count - 1) * 6)

    [void]$svg.AppendLine("<line x1='$($plotX - 6)' y1='$y' x2='$plotX' y2='$y' stroke='#000' stroke-width='0.8'/>")
    [void]$svg.AppendLine("<text class='term' x='$labelX' y='$labelY' text-anchor='end'>")
    for ($lineIndex = 0; $lineIndex -lt $wrappedLines.Count; $lineIndex++) {
      $dy = if ($lineIndex -eq 0) { 0 } else { 12 }
      [void]$svg.AppendLine("<tspan x='$labelX' dy='$dy'>$(Escape-Xml $wrappedLines[$lineIndex])</tspan>")
    }
    [void]$svg.AppendLine("</text>")
    [void]$svg.AppendLine("<circle cx='$x' cy='$y' r='$radius' fill='$color' stroke='none'/>")
  }

  $legendX = $plotX + $plotWidth + 42
  $legendY = $plotY + 52
  [void]$svg.AppendLine("<text class='legend' x='$legendX' y='$($legendY - 16)' font-weight='700'>p.adjust</text>")
  [void]$svg.AppendLine("<rect x='$legendX' y='$legendY' width='17' height='132' fill='url(#plasma_adjust_gradient)'/>")

  [void]$svg.AppendLine("<text class='legend' x='$($legendX + 24)' y='$($legendY + 5)'>$([Math]::Round($minFdr, 3))</text>")
  [void]$svg.AppendLine("<text class='legend' x='$($legendX + 24)' y='$($legendY + 131)'>$([Math]::Round($maxFdr, 3))</text>")

  [void]$svg.AppendLine("<text class='legend' x='$legendX' y='$($legendY + 178)' font-weight='700'>Count</text>")
  $countValues = @($minCount, [Math]::Round(($minCount + $maxCount) / 2), $maxCount) | Select-Object -Unique

  for ($countIndex = 0; $countIndex -lt $countValues.Count; $countIndex++) {
    $countValue = [double]$countValues[$countIndex]
    $countScale = if ($maxCount -le $minCount) { 0.5 } else { ($countValue - $minCount) / ($maxCount - $minCount) }
    $radius = 3 + (9 * [Math]::Sqrt([Math]::Max(0, $countScale)))
    $circleY = $legendY + 204 + ($countIndex * 35)
    [void]$svg.AppendLine("<circle cx='$($legendX + 9)' cy='$circleY' r='$radius' fill='#333'/>")
    [void]$svg.AppendLine("<text class='legend' x='$($legendX + 28)' y='$($circleY + 4)'>$([int]$countValue)</text>")
  }

  [void]$svg.AppendLine("</g>")
}

[void]$svg.AppendLine("</svg>")
Set-Content -Path $outputSvg -Value $svg.ToString() -Encoding UTF8
Write-Output "Saved SVG: $outputSvg"
