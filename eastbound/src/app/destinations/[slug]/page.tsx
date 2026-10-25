import type { Metadata } from "next";
import Image from "next/image";
import Link from "next/link";
import { notFound } from "next/navigation";
import { destinations, getDestination } from "@/content/destinations";

type Props = { params: Promise<{ slug: string }> };

export function generateStaticParams() {
  return destinations.map((d) => ({ slug: d.slug }));
}

export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const { slug } = await params;
  const d = getDestination(slug);
  if (!d) return { title: "Destination" };
  return {
    title: d.name,
    description: d.blurb,
  };
}

export default async function DestinationPage({ params }: Props) {
  const { slug } = await params;
  const d = getDestination(slug);
  if (!d) notFound();

  return (
    <>
      <div className="cover-band">
        <Image src={d.image} alt={d.imageAlt} fill priority sizes="100vw" />
        <div className="cover-band__shade" />
        <div className="cover-band__label wrap" style={{ left: 0, right: 0 }}>
          <p className="eyebrow">
            {d.region} · {d.days}
          </p>
          <h1 className="prose-title" style={{ marginBottom: 0 }}>
            {d.name}
          </h1>
        </div>
      </div>

      <div className="wrap section--tight" style={{ paddingTop: 0 }}>
        <p className="muted" style={{ maxWidth: "40rem", fontSize: "1.1rem" }}>
          {d.blurb}
        </p>
        <p className="meta-row">
          <span>Best for: {d.bestFor}</span>
          <span>Suggested stay: {d.days}</span>
        </p>

        <div className="split">
          <section className="stack-gaps">
            <div>
              <h2 style={{ fontFamily: "var(--font-display)", fontWeight: 560 }}>See</h2>
              <ul className="plain-list">
                {d.see.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            </div>
            <div>
              <h2 style={{ fontFamily: "var(--font-display)", fontWeight: 560 }}>Eat</h2>
              <ul className="plain-list">
                {d.eat.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            </div>
          </section>
          <aside className="stack-gaps">
            <div className="panel">
              <h3>Move</h3>
              <ul className="plain-list">
                {d.move.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            </div>
            <div className="panel">
              <h3>Stay</h3>
              <p className="muted" style={{ marginBottom: 0 }}>
                {d.stay}
              </p>
            </div>
            <Link href="/transport" className="btn btn--ghost">
              Transport how-to
            </Link>
          </aside>
        </div>

        <p style={{ marginTop: "2.5rem" }}>
          <Link href="/destinations">← All destinations</Link>
        </p>
      </div>
    </>
  );
}
