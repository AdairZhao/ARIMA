import type { Metadata } from "next";
import Link from "next/link";

export const metadata: Metadata = {
  title: "Essentials",
  description:
    "Alipay, WeChat Pay, eSIM, maps, and etiquette for foreign visitors to China.",
};

const blocks = [
  {
    title: "Payments",
    body: "Set up Alipay and WeChat Pay with an eligible foreign card before departure. Test a small payment at home if possible. Keep a little cash for rare offline stalls, but expect QR everywhere.",
  },
  {
    title: "Connectivity",
    body: "A China-capable eSIM (or local SIM) is non-negotiable for maps and ride-hail. Hotel Wi‑Fi alone will strand you. Download offline translation packs.",
  },
  {
    title: "Maps that work",
    body: "Apple Maps and Amap (Gaode) are more dependable on the ground than Google Maps. Save hotel pins and metro exits while online.",
  },
  {
    title: "Tickets & passport",
    body: "Major museums and parks often need advance real-name booking. Carry your passport daily — hotels, trains, and attractions will ask.",
  },
  {
    title: "Etiquette quick hits",
    body: "Queues matter at stations. Speak softly on trains. Photos of people: ask. Tap water is not for drinking — hotels provide kettles and sealed bottles.",
  },
];

export default function EssentialsPage() {
  return (
    <div className="wrap section--tight">
      <header className="page-hero">
        <p className="eyebrow">Essentials</p>
        <h1 className="prose-title">Survival kit for day one</h1>
        <p>
          If you only prepare five things, make them payments, data, maps,
          passport tickets, and a calm plan for arrival night.
        </p>
      </header>

      <div className="check-list">
        {blocks.map((b) => (
          <article key={b.title} className="check-item">
            <h2 style={{ fontFamily: "var(--font-display)", fontSize: "1.45rem", margin: 0 }}>
              {b.title}
            </h2>
            <p className="muted" style={{ margin: 0 }}>
              {b.body}
            </p>
          </article>
        ))}
      </div>

      <div className="cta-band" style={{ marginTop: "3rem" }}>
        <h2>Put it on a timeline</h2>
        <p>The plan checklist orders these tasks by when to do them.</p>
        <Link href="/plan" className="btn btn--primary">
          Go to plan checklist
        </Link>
      </div>
    </div>
  );
}
