import pytest
from pubmed_search import PubmedSearch

class TestPub():
    title = 'Arabidopsis synchronizes jasmonate-mediated defense with insect circadian behavior.'
    authors = 'Goodspeed D, Chehab EW, Min-Venditti A, Braam J, Covington MF'
    journal = 'Proc Natl Acad Sci U S A'
    year = '2012'
    month = 3
    day = '20'
    url = 'https://www.pnas.org/content/109/12/4674'
    pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed/22331878'
    citation = 'Goodspeed D, Chehab EW, Min-Venditti A, Braam J, Covington MF (2012). Arabidopsis synchronizes jasmonate-mediated defense with insect circadian behavior. Proc Natl Acad Sci U S A 109(12): 4674-7.'
    minicitation = 'Goodspeed D - Covington MF - 2012 - Proc Natl Acad Sci U S A'
    abstract = "Diverse life forms have evolved internal clocks enabling them to monitor time and thereby anticipate the daily environmental changes caused by Earth's rotation. The plant circadian clock regulates expression of about one-third of the Arabidopsis genome, yet the physiological relevance of this regulation is not fully understood. Here we show that the circadian clock, acting with hormone signals, provides selective advantage to plants through anticipation of and enhanced defense against herbivory. We found that cabbage loopers (Trichoplusia ni) display rhythmic feeding behavior that is sustained under constant conditions, and plants entrained in light/dark cycles coincident with the entrainment of the T. ni suffer only moderate tissue loss due to herbivory. In contrast, plants entrained out-of-phase relative to the insects are significantly more susceptible to attack. The in-phase entrainment advantage is lost in plants with arrhythmic clocks or deficient in jasmonate hormone; thus, both the circadian clock and jasmonates are required. Circadian jasmonate accumulation occurs in a phase pattern consistent with preparation for the onset of peak circadian insect feeding behavior, providing evidence for the underlying mechanism of clock-enhanced herbivory resistance. Furthermore, we find that salicylate, a hormone involved in biotrophic defense that often acts antagonistically to jasmonates, accumulates in opposite phase to jasmonates. Our results demonstrate that the plant circadian clock provides a strong physiological advantage by performing a critical role in Arabidopsis defense."

@pytest.fixture
def example_pubmed_search():
    return PubmedSearch('gixos30110@lowdh.com') #fake temp mail

@pytest.mark.parametrize("search_query, start, max", [
    ('cancer', 0, 1),
    ('autism', 2, 18),
])

def test_get_ids(example_pubmed_search, search_query, start, max):
    result = example_pubmed_search.get_ids(search_query, start, max)
    assert type(result) is list
    assert result[0] == example_pubmed_search.get_ids(search_query, 0, start+1)[-1]
    assert len(result) == max

def test_get_publication(example_pubmed_search):
    test_pub = TestPub()
    pub = example_pubmed_search.get_publication('22331878')
    assert pub.title == test_pub.title
    assert pub.authors == test_pub.authors
    assert pub.journal == test_pub.journal
    assert pub.year == test_pub.year
    assert pub.month == test_pub.month
    assert pub.day == test_pub.day
    assert pub.url == test_pub.url
    assert pub.pubmed_url == test_pub.pubmed_url
    assert pub.cite() == test_pub.citation
    assert pub.cite_mini() == test_pub.minicitation
    assert pub.abstract == test_pub.abstract