import type { Metadata } from "next";
import Image from "next/image";
import Link from "next/link";
import { notFound } from "next/navigation";
import { getItinerary, itineraries } from "@/content/itineraries";

type Props = { params: Promise<{ slug: string }> };

export function generateStaticParams() {
  return itineraries.map((i) => ({ slug: i.slug }));
}

export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const { slug } = await params;
  const item = getItinerary(slug);
  if (!item) return { title: "Itinerary" };
  return { title: item.title, description: item.summary };
}

export default async function ItineraryPage({ params }: Props) {
  const { slug } = await params;
  const item = getItinerary(slug);
  if (!item) notFound();

  return (
    <>
      <div className="cover-band">
        <Image src={item.image} alt={item.imageAlt} fill priority sizes="100vw" />
        <div className="cover-band__shade" />
        <div className="cover-band__label wrap" style={{ left: 0, right: 0 }}>
          <p className="eyebrow">{item.days}-day template</p>
          <h1 className="prose-title" style={{ marginBottom: 0 }}>
            {item.title}
          </h1>
        </div>
      </div>

      <div className="wrap section--tight" style={{ paddingTop: 0 }}>
        <p className="muted" style={{ maxWidth: "40rem", fontSize: "1.1rem" }}>
          {item.summary}
        </p>
        <div className="route-chips">
          {item.route.map((city) => (
            <span key={city}>{city}</span>
          ))}
        </div>

        <section>
          {item.daysPlan.map((block) => (
            <article key={block.day} className="day-block">
              <h2 style={{ fontFamily: "var(--font-display)", fontSize: "1.4rem", margin: 0 }}>
                {block.day} · {block.focus}
              </h2>
              <p className="muted" style={{ margin: "0.4rem 0 0" }}>
                {block.notes}
              </p>
            </article>
          ))}
        </section>

        <div className="cta-band" style={{ marginTop: "2.5rem" }}>
          <h2>Pair this with the pre-trip checklist</h2>
          <p>Visa, payments, and apps before you lock hotel nights.</p>
          <Link href="/plan" className="btn btn--primary">
            Open plan checklist
          </Link>
        </div>
      </div>
    </>
  );
}
