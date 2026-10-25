import type { Metadata } from "next";
import Image from "next/image";
import Link from "next/link";
import { itineraries } from "@/content/itineraries";

export const metadata: Metadata = {
  title: "Itineraries",
  description: "Ready-made China trip templates for first-time visitors.",
};

export default function ItinerariesPage() {
  return (
    <div className="wrap section--tight">
      <header className="page-hero">
        <p className="eyebrow">Itineraries</p>
        <h1 className="prose-title">Routes that already fit together</h1>
        <p>
          Fixed lengths, realistic train hops, and days that leave room for jet
          lag — not fantasy packing.
        </p>
      </header>

      <div className="destination-rail" style={{ gridTemplateColumns: "1fr" }}>
        {itineraries.map((item) => (
          <Link
            key={item.slug}
            href={`/itineraries/${item.slug}`}
            className="dest-link"
            style={{ minHeight: 320 }}
          >
            <Image src={item.image} alt={item.imageAlt} fill sizes="100vw" />
            <div className="dest-link__body">
              <p className="dest-link__meta">{item.days} days</p>
              <h3>{item.title}</h3>
              <p className="dest-link__meta">{item.summary}</p>
            </div>
          </Link>
        ))}
      </div>
    </div>
  );
}
