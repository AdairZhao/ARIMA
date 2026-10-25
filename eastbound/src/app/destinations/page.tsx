import type { Metadata } from "next";
import Image from "next/image";
import Link from "next/link";
import { destinations } from "@/content/destinations";

export const metadata: Metadata = {
  title: "Destinations",
  description: "English city guides for Beijing, Shanghai, Xi'an, and Chengdu.",
};

export default function DestinationsPage() {
  return (
    <div className="wrap section--tight">
      <header className="page-hero">
        <p className="eyebrow">Destinations</p>
        <h1 className="prose-title">Hub cities that teach you China</h1>
        <p>
          Begin with places that have trains, tickets, and enough English
          support — then branch out.
        </p>
      </header>

      <div className="destination-rail">
        {destinations.map((d) => (
          <Link key={d.slug} href={`/destinations/${d.slug}`} className="dest-link">
            <Image src={d.image} alt={d.imageAlt} fill sizes="(max-width: 720px) 100vw, 25vw" />
            <div className="dest-link__body">
              <p className="dest-link__meta">
                {d.region} · {d.days}
              </p>
              <h3>{d.name}</h3>
              <p className="dest-link__meta">{d.blurb}</p>
            </div>
          </Link>
        ))}
      </div>
    </div>
  );
}
